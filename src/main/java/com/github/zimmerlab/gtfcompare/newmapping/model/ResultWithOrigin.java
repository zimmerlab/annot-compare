package com.github.zimmerlab.gtfcompare.newmapping.model;

import com.github.zimmerlab.gtfcompare.model.GenePair;

public record ResultWithOrigin(GenePair genePair, MappingOrigin origin) {
}
