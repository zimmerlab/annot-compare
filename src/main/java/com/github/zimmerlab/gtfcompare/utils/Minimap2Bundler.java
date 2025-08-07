package com.github.zimmerlab.gtfcompare.utils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

public class Minimap2Bundler {
    public static Path extractMinimap2() throws IOException {
        var os = System.getProperty("os.name").toLowerCase();
        var arch = System.getProperty("os.arch").toLowerCase();
        String platform;
        if (os.contains("linux") && arch.contains("64")) {
            platform = "linux-amd64";
        } else if (os.contains("mac")) {
            platform = "macos-amd64";
            throw new UnsupportedOperationException("Unsupported OS/Arch: " + os + "/" + arch);
        } else if (os.contains("win")) {
            platform = "win-amd64";
            throw new UnsupportedOperationException("Unsupported OS/Arch: " + os + "/" + arch);
        } else {
            throw new UnsupportedOperationException("Unsupported OS/Arch: " + os + "/" + arch);
        }
        var resourcePath = "/minimap2/" + platform + "/minimap2" + (platform.startsWith("win") ? ".exe" : "");

        try (InputStream in = Minimap2Bundler.class.getResourceAsStream(resourcePath)) {
            if (in == null) {
                throw new FileNotFoundException("Resource not found: " + resourcePath);
            }
            var exe = Files.createTempFile("minimap2-", platform.startsWith("win") ? ".exe" : "");
            exe.toFile().deleteOnExit();
            Files.copy(in, exe, StandardCopyOption.REPLACE_EXISTING);
            exe.toFile().setExecutable(true, true);
            return exe;
        }
    }
}
