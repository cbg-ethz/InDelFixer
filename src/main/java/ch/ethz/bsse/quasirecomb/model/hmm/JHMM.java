package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorker;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorkerRecalc;
import ch.ethz.bsse.quasirecomb.utils.Random;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.HashMap;
import java.util.Map;

/**
 * Offers a HMM with the E and M step.
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class JHMM {

    private final int N;
    private final int L;
    private final int K;
    private final int n;
    //rho[j][l][k] := transition prob. at position j, for k given l
    private double[][][] rho;
    private double[][][] priorRho;
    private double[] pi;
    private double[][][] mu;
    private double eps;
    private double antieps;
    private double[][] nJK;
    private double[][][] nJKL;
    private double[][][] nJKV;
    private double[][] nVB;
    private double[] nJeq;
    private double[] nJneq;
    private double neq;
    private double nneq;
    private double loglikelihood;
    private double likelihood;
    private double Q;
    private Map<ReadHMM, Integer> readHMMMap;
    private int mapSize;
    private Map<byte[], Integer> clusterReads;

    public JHMM(String[] reads, int K, int n) {
        this(Utils.splitReadsIntoByteArrays(reads), K, n, Globals.ESTIMATION_EPSILON);
    }

    public JHMM(String[] reads, int K, int n, double epsilon) {
        this(Utils.splitReadsIntoByteArrays(reads), K, n, epsilon);
    }

    public JHMM(byte[][] reads, int K, int n) {
        this(reads, K, n, Globals.ESTIMATION_EPSILON);
    }

    public JHMM(byte[][] reads, int K, int n, double epsilon) {
        this(reads, reads.length, reads[0].length, K, n, epsilon);
    }

    public JHMM(Map<byte[], Integer> reads, int N, int L, int K, int n, OptimalResult or) {
        this(reads, N, L, K, n, Globals.ESTIMATION_EPSILON, (1 - (n - 1) * Globals.ESTIMATION_EPSILON), or.getRho(), or.getPi(), or.getMu(), or.getPriorRho());
    }

    public JHMM(Map<byte[], Integer> reads, int N, int L, int K, int n, double epsilon) {
        this(reads, N, L, K, n, epsilon,
                (1 - (n - 1) * epsilon),
                Random.generateInitRho(L - 1, K),
                Random.generateInitPi(K),
                Random.generateMuInit(L, K, n),
                Random.generatePriorRho(L - 1, K));
    }

    public JHMM(byte[][] reads, int N, int L, int K, int n, double epsilon) {
        this(reads, N, L, K, n, epsilon,
                (1 - (n - 1) * epsilon),
                Random.generateInitRho(L - 1, K),
                Random.generateInitPi(K),
                Random.generateMuInit(L, K, n),
                Random.generatePriorRho(L - 1, K));
    }

    public JHMM(byte[][] reads, int N, int L, int K, int n, double eps, double antieps, double[][][] rho, double[] pi, double[][][] mu, double[][][] priorRho) {
        this.N = N;
        this.L = L;
        this.K = K;
        this.n = n;
        this.clusterReads = Utils.clusterReads(reads);
        this.rho = rho;
        this.mu = mu;
        this.eps = eps;
        this.antieps = antieps;
        this.pi = pi;
        this.priorRho = priorRho;
        this.start();
        this.calculate();
    }

    private JHMM(Map<byte[], Integer> reads, int N, int L, int K, int n, double eps, double antieps, double[][][] rho, double[] pi, double[][][] mu, double[][][] priorRho) {
        this.N = N;
        this.L = L;
        this.K = K;
        this.n = n;
        this.clusterReads = reads;
        this.rho = rho;
        this.mu = mu;
        this.eps = eps;
        this.antieps = antieps;
        this.pi = pi;
        this.priorRho = priorRho;
        this.start();
        this.calculate();
    }

    private void calculateLoglikelihood() {
        long time = System.currentTimeMillis();
        this.loglikelihood = 0d;
        for (ReadHMM r : this.readHMMMap.keySet()) {
            for (int j = 0; j < L; j++) {
                this.loglikelihood += Math.log(r.getC(j)) * this.readHMMMap.get(r);
            }
        }
        this.likelihood = 0d;
        for (ReadHMM r : this.readHMMMap.keySet()) {
            for (int i = 0; i < this.readHMMMap.get(r); i++) {
                double llh = 1d;
                for (int j = 0; j < L; j++) {
                    llh *= r.getC(j);
                }
                this.likelihood += llh;
            }
        }
        System.out.println("L\t: " + (System.currentTimeMillis() - time));
    }

    private void start() {
        long time = System.currentTimeMillis();
        byte[][] uniqueReads = new byte[clusterReads.keySet().size()][L];
        for (byte[] read : this.clusterReads.keySet()) {
            uniqueReads[mapSize++] = read;
        }
        this.readHMMMap = new HashMap<>();
        ReadHMM[] rArray = new ReadHMM[uniqueReads.length];
        if (Globals.PARALLEL_JHMM) {
            rArray = Globals.fjPool.invoke(new ReadHMMWorker(this, uniqueReads, rho, pi, mu, eps, antieps, K, L, n, 0, mapSize)).toArray(new ReadHMM[mapSize]);
        } else {
            for (int i = 0; i < uniqueReads.length; i++) {
                rArray[i] = new ReadHMM(L, K, n, uniqueReads[i], rho, pi, mu, eps, antieps);
            }
        }
        for (ReadHMM r : rArray) {
            this.readHMMMap.put(r, clusterReads.get(r.getRead()));
        }
        System.out.println("R\t: " + (System.currentTimeMillis() - time));
    }

    public void restart() {
        long time = System.currentTimeMillis();
        if (Globals.PARALLEL_JHMM) {
            Globals.fjPool.invoke(new ReadHMMWorkerRecalc(this.readHMMMap.keySet().toArray(new ReadHMM[this.readHMMMap.keySet().size()]), rho, pi, mu, eps, 0, mapSize));
        } else {
            for (ReadHMM r : this.readHMMMap.keySet()) {
                r.recalc(rho, pi, mu, eps);
            }
        }
        System.out.println("R\t: " + (System.currentTimeMillis() - time));
        this.calculate();
    }

    private void calculate() {
        this.eStep();
        this.calculateLoglikelihood();
        this.mStep();
    }

    private void eStep() {
        long time = System.currentTimeMillis();
        this.nJK = new double[L][K];
        this.nJKL = new double[L][K][K];
        this.nJKV = new double[L][K][n];
        this.nVB = new double[n][n];
        this.nJeq = new double[L];
        this.nJneq = new double[L];
        for (ReadHMM r : this.readHMMMap.keySet()) {
            for (int x = 0; x < this.readHMMMap.get(r); x++) {
                for (int j = 0; j < L; j++) {
                    for (int k = 0; k < K; k++) {
                        this.nJK[j][k] += r.gamma(j, k);
                        if (j > 0) {
                            for (int l = 0; l < K; l++) {
                                this.nJKL[j][k][l] += r.xi(j, k, l);
                                if (k == l) {
                                    this.nJeq[j] += r.xi(j, k, l);
                                } else {
                                    this.nJneq[j] += r.xi(j, k, l);
                                }
                            }
                        }
                        for (int v = 0; v < n; v++) {
                            this.nJKV[j][k][v] += r.gamma(j, k, v);
                        }
                    }
                    for (int v = 0; v < n; v++) {
                        for (int b = 0; b < n; b++) {
                            if (r.getRead()[j] == b) {
                                for (int k = 0; k < K; k++) {
                                    this.nVB[v][b] += r.gamma(j, k, v);
                                    if (v == b) {
                                        this.neq += r.gamma(j, k, v);
                                    } else {
                                        this.nneq += r.gamma(j, k, v);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.println("E\t: " + (System.currentTimeMillis() - time));
    }

    private double[][][] calcMu() {
        double[][][] mu_ = new double[L][K][n];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double sumV = 0d;
                for (int v = 0; v < n; v++) {
                    mu_[j][k][v] = this.rho_f(this.getnJKV(j, k, v) + Globals.PRIOR_ALPHA);
                    sumV += mu_[j][k][v];
                }
                if (sumV != 0) {
                    for (int v = 0; v < n; v++) {
                        mu_[j][k][v] /= sumV;
                    }
                } else {
                    System.out.println("mu");
                }
            }
        }
        return mu_;
    }

    private double[][][] calcRho() {
        double[][][] rho_ = new double[L - 1][K][K];
        for (int j = 1; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double divisor = 0d;
                for (int l = 0; l < K; l++) {
                    if (Double.isNaN(this.getnJKL(j, k, l))) {
                        System.out.println("rho");
                    }
                    rho_[j - 1][k][l] = this.rho_f(this.getnJKL(j, k, l) + Globals.PRIOR_ALPHA);
                    divisor += rho_[j - 1][k][l];
                }

//                for (int l = 0; l < K; l++) {
//                    divisor += this.getnJKL(j, k, l);
//                }
//                divisor = this.rho_f(divisor + Globals.PRIOR_ALPHA);
                for (int l = 0; l < K; l++) {
                    rho_[j - 1][k][l] /= divisor;
                }
            }
        }
        return rho_;
    }

    private double rho_f(double upsilon) {
        if (upsilon == 0d) {
            return 0d;
        }
//        System.out.println("#"+upsilon);
        return Math.exp(rho_phi(upsilon));
    }

    private double rho_phi(double upsilon) {
        double x = -1d;
        try { 
            x = (upsilon > 7) ? rho_g(upsilon - .5) : (rho_phi(upsilon + 1) - 1 / upsilon);
            
        } catch (StackOverflowError e) {
            System.err.println(upsilon);
            System.err.println(e);
        }
        return x;
    }

    private double rho_g(double x) {
        return Math.log(x) + .04167 * Math.pow(x, -2) - .00729 * Math.pow(x, -4) + .00384 * Math.pow(x, -6) - .00413 * Math.pow(x, -8);
    }

    private double[] calcPi() {
        double[] pi_ = new double[K];
        double sumK = 0d;
        for (int k = 0; k < K; k++) {
            pi_[k] = this.getnJK(0, k);
            sumK += pi_[k];
        }
        for (int k = 0; k < K; k++) {
            pi_[k] /= sumK;
        }
        return pi_;
    }

    private void mStep() {
        long time = System.currentTimeMillis();
        if (!Globals.rho0) {
            this.rho = this.calcRho();
        }
        this.pi = this.calcPi();
        this.mu_old = this.mu;
        this.mu = this.calcMu();
        System.out.println("M\t: " + (System.currentTimeMillis() - time));
//        this.eps = (1d/(n-1d))*(this.nneq/((double)this.nneq+this.neq));
        System.out.println("#EPS: "+eps);
    }
    public double[][][] mu_old;

    public double getnJK(int j, int k) {
        return nJK[j][k];
    }

    public double getnJKL(int j, int k, int l) {
        return nJKL[j][k][l];
    }

    public double getnJKV(int j, int k, int v) {
        return nJKV[j][k][v];
    }

    public double getnJeq(int j) {
        return nJeq[j];
    }

    public double getnJneq(int j) {
        return nJneq[j];
    }

    public double[][] getnVB() {
        return nVB;
    }

    public double getNeq() {
        return neq;
    }

    public double getNneq() {
        return nneq;
    }

    public int getK() {
        return K;
    }

    public int getL() {
        return L;
    }

    public int getN() {
        return N;
    }

    public double getAntieps() {
        return antieps;
    }

    public double getEps() {
        return eps;
    }

    public int getn() {
        return n;
    }

    public double getLoglikelihood() {
        return loglikelihood;
    }

    public double getQ() {
        return Q;
    }

    public double[][][] getMu() {
        return mu;
    }

    public double[] getPi() {
        return pi;
    }

    public double[][][] getRho() {
        return rho;
    }

    public Map<ReadHMM, Integer> getReadHMMMap() {
        return readHMMMap;
    }

    public double getLikelihood() {
        return likelihood;
    }

    public double[][][] getPrior_rho() {
        return priorRho;
    }
}
