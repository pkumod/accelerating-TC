We use five state-of-the-art implementations in our paper, and we implement four by ourselves, including Hu's work, Bission's work, Fox's work, Tricore. Those codes could be make by Makefile in them. And we use gunrock's published version in https://github.com/gunrock/gunrock.

The source code of our A-direction strategy is in folder "edge_direction". We implement our algorithm in both serilized and parallelized manner, and the code can be compiled directly using g++ and CUDA compiler nvcc.

The source code of our A-order strategy and four other state-of-the-art reordering methods are in folder "vertex_ordering".

The folder "tools" contains some data preprocessing codes.


