package com.github.zimmerlab.gtfcompare.utils;

import org.springframework.util.StopWatch;

import java.util.ArrayList;
import java.util.List;

public class Constants {
    public final static List<StopWatch> STOP_WATCHES = new ArrayList<>();

    // Comparator Names
    public final static String LENGTH_COMPARATOR_NAME = "length";
    public final static String SEQUENCE_COMPARATOR_NAME = "sequence";
    public final static String START_COMPARATOR_NAME = "start";
    public final static String STOP_COMPARATOR_NAME = "stop";
}
