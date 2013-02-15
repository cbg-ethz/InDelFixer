///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package ch.ethz.bsse.indelfixer.refinement;
//
//import ch.ethz.bsse.indelfixer.stored.Read;
//import com.google.common.collect.Lists;
//import java.io.File;
//import java.util.HashMap;
//import java.util.List;
//import java.util.Map;
//import java.util.concurrent.ArrayBlockingQueue;
//import java.util.concurrent.BlockingQueue;
//import java.util.concurrent.ExecutionException;
//import java.util.concurrent.ExecutorService;
//import java.util.concurrent.Future;
//import java.util.concurrent.RejectedExecutionHandler;
//import java.util.concurrent.ThreadPoolExecutor;
//import java.util.concurrent.TimeUnit;
//import net.sf.samtools.SAMFileReader;
//import net.sf.samtools.SAMRecord;
//
///**
// *
// * @author toepfera
// */
//public class Refining {
//
//    private static final BlockingQueue<Runnable> blockingQueue = new ArrayBlockingQueue<>(Runtime.getRuntime().availableProcessors() - 1);
//    private static final RejectedExecutionHandler rejectedExecutionHandler = new ThreadPoolExecutor.CallerRunsPolicy();
//    private static ExecutorService executor = refreshExecutor();
//
//    private static ExecutorService refreshExecutor() {
////        return Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() - 1);
//        return new ThreadPoolExecutor(Runtime.getRuntime().availableProcessors() - 1, Runtime.getRuntime().availableProcessors() - 1, 0L, TimeUnit.MILLISECONDS, blockingQueue, rejectedExecutionHandler);
//    }
//
//    public static Read[] parseBAMSAM(String location) {
//        Map<String, ReadTMP> readMap = new HashMap<>();
//        File bam = new File(location);
//        SAMFileReader sfr = new SAMFileReader(bam);
//        double size = 0;
//        for (final SAMRecord samRecord : sfr) {
//            size++;
//        }
//        sfr.close();
//        sfr = new SAMFileReader(bam);
//        List<Future<ReadTMP>> readFutures = Lists.newArrayListWithExpectedSize((int) size);
//        int counter = 0;
//        for (final SAMRecord samRecord : sfr) {
//            readFutures.add(executor.submit(new SFRComputing(samRecord)));
//            Globals.getINSTANCE().print("Parsing\t\t" + (Math.round((counter++ / size) * 100)) + "%");
//        }
//        Globals.getINSTANCE().print("Parsing\t\t100%");
//        sfr.close();
//        for (Future<ReadTMP> future : readFutures) {
//            try {
//                ReadTMP read = future.get();
//                if (read != null) {
//                    String name = read.name;
//                    int refStart = read.refStart;
//                    byte[] readBases = read.readBases;
//                    double[] quality = read.quality;
//                    boolean hasQuality = read.hasQuality;
//                    if (readMap.containsKey(name)) {
//                        if (hasQuality) {
//                            readMap.get(name).setPairedEnd(BitMagic.pack(readBases), refStart, refStart + readBases.length, quality);
//                        } else {
//                            readMap.get(name).setPairedEnd(BitMagic.pack(readBases), refStart, refStart + readBases.length);
//                        }
//                        Read r2 = readMap.get(name);
//                        if (r2.getCrickBegin() - r2.getWatsonEnd() > 250) {
//                            System.out.println(name + "\t" + (r2.getCrickBegin() - r2.getWatsonEnd()));
//                            readMap.put(name, null);
//                        }
//                        if (Globals.getINSTANCE().isUNPAIRED()) {
//                            Read r = readMap.get(name);
//                            if (r.isPaired()) {
//                                r.unpair();
//                                if (hasQuality) {
//                                    readMap.put(name + "_R", new Read(BitMagic.pack(readBases), refStart, refStart + readBases.length, quality));
//                                } else {
//                                    readMap.put(name + "_R", new Read(BitMagic.pack(readBases), refStart, refStart + readBases.length));
//                                }
//                            }
//                        }
//                        Globals.getINSTANCE().incPAIRED();
//                    } else {
//                        if (hasQuality) {
//                            readMap.put(name, new Read(BitMagic.pack(readBases), refStart, refStart + readBases.length, quality));
//                        } else {
//                            readMap.put(name, new Read(BitMagic.pack(readBases), refStart, refStart + readBases.length));
//                        }
//                    }
//                }
//            } catch (InterruptedException | ExecutionException ex) {
//                System.err.println(ex);
//            }
//        }
//        readFutures.clear();
//
//        Map<Integer, Read> hashed = new HashMap<>();
//        for (Read r1 : readMap.values()) {
//            if (r1 != null) {
//                int hash = r1.hashCode();
//                if (hashed.containsKey(hash)) {
//                    hashed.get(hash).incCount();
//                } else {
//                    hashed.put(hash, r1);
//                }
//            }
//        }
////        return null;
//        return hashed.values().toArray(new Read[hashed.size()]);
//    }
//}
