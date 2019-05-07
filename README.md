# Openmp_Google_Page_Rank_Team-11
 This project is an openmp implementation of Google Page Rank algorithm.Our approach is an implementation of the research paper whose details are mentioned in the project report.
 
 The algorithm used follows two kind of comparisons:
 
 a) Baseline Algo(Iterative) vs Parallel Baseline Algo: The baseline algorithm is an iterative approach to find the pagerank ,while the parallel Baseline Algo is just the openmp implementation of the baseline algorithm.
    To run the baseline algorithm following format is to be used:
         
         1.To compile use the following format:
            g++ parallel_pgrank_baseline.cpp -fopenmp
            
         2.To run:     
            ./a.out TestFilename
            eg: ./a.out < web-Stanford.txt


b) Optimized Algo(Iterative) vs Parallel Optimized Algo: In optimized algorithm we first find strongly connected components in the graph and then use topological sort to divide graph into different levels as described in report.Now in parallel version we can run all strongly connected components within each level at same time.
        
        1. To compile use the following format:
            g++ parallel_pgrank_optimized.cpp -fopenmp

        2. To run:
           ./a.out TestFilename
           eg: ./a.out < web-Stanford.txt

