This is an extension of the id3 program I wrote in decision_tree. This has two parts. First, improve the id3 algorithm I wrote. Second, add the -t option to the program.

If you call the program without the -t it will perform as it did in decision_tree. With one exception. When about to recurse on the ID3Build with no data instead return an Ans using maxAns for the current data rather than the subsetted data which is empty for that value of the feature. At the end of printing all the output for the tree as before have it also print the actual tree data structure with print(tree).

If the -t option is supplied have it read in the tree data structure from the file whose name is given with the -t option. Code for reading in the file and putting the data structure will be supplied. Now proceed to read from standard input as before using readProblem(). For each read data vector print the data including the expected answer and then the answer the tree would give.
