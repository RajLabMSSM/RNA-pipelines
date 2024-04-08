# Tiebrush pipeline

Runs tiebrush and tiecov. Shrinks down the size of BAM files by collapsing duplicate reads.

Tiecov creates BigWig and junction files for viewing coverage and splice junctions respectively.

To run, create a set of folders, one folder per merged group. Then symlink the BAM and BAI files all the samples in that group

Tiebrush will then aggregate those BAMs together to produce a grouped BAM, one per directory.
