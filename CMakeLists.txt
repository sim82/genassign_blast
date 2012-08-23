include_directories( include )

add_library( blib src/alignments.c src/blast.c src/bytepackGappedScoring.c src/descriptions.c src/encoding.c src/fasterGappedExtension.c src/gappedExtension.c src/gappedScoring.c src/nuGappedScoring.c src/global.c src/hitMatrix.c src/karlin.c src/memBlocks.c src/memSingleBlock.c src/nucleotideLookup.c src/oldGappedScoring.c src/oldSemiGappedScoring.c src/parameters.c src/print.c src/PSSMatrix.c src/qPosList.c src/readFasta.c src/readFile.c src/scoreMatrix.c src/semiGappedScoring.c src/statistics.c src/ungappedExtension.c src/wordLookupDFA.c src/writeFile.c src/constants.c src/smithWatermanTraceback.c src/smithWatermanScoring.c src/tableGappedScoring.c src/vbyte.c src/unpack.c src/index.c src/postings.c src/hashcounter.c src/writedb.c src/readdb.c src/search.c src/wildcards.c src/dust.c src/seg.c )

add_executable( blast src/blast.c )
target_link_libraries( blast blib m )

add_executable( formatdb src/formatdb.c )
target_link_libraries( formatdb blib m )

