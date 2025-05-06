package com.github.zimmerlab.gtfcompare.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.HashMap;

public class ProgramEvaluator {
    private final static Logger LOG  = LogManager.getLogger(ProgramEvaluator.class);
    public static void evaluate(String out){
        var taskTimeMap = new HashMap<String, Long>();

        for(var stopWatch : Constants.STOP_WATCHES){
            for(var task : stopWatch.getTaskInfo()){
                var taskName = task.getTaskName();
                var value = taskTimeMap.getOrDefault(taskName, 0L);
                taskTimeMap.put(taskName, value + task.getTimeMillis());
            }
        }

        var now = LocalDateTime.now();
        var formatter = DateTimeFormatter.ofPattern("yyyyMMddHHmmss");
        var formattedDate = now.format(formatter);
        String fileName ;
        fileName = "time_eval_" + formattedDate + ".tsv";


        try(var writer = new BufferedWriter(new FileWriter(new File(out, fileName)))){
            writer.write("Task Name\tMilliseconds\n");
            for(var task : taskTimeMap.entrySet()){
                System.out.println(task.getKey() + ": " + task.getValue() + " ms");
                writer.write(task.getKey() + "\t" + task.getValue() + "\n");
            }
        } catch (Exception e){
            LOG.error("Error while printing stats", e);
        }

    }
}
