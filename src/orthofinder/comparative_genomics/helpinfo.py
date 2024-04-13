from ..utils import util
from othofinder import nThreadsDefault


def PrintHelp():
    print("Usage")    
    print("-----")
    print("orthologues.py orthofinder_results_directory [-t max_number_of_threads]")
    print("orthologues.py -h")
    print("\n")
    
    print("Arguments")
    print("---------")
    print("""orthofinder_results_directory
    Generate gene trees for the orthogroups, generated rooted species tree and infer ortholgues.\n""")
    
    print(("""-t max_number_of_threads, --threads max_number_of_threads
    The maximum number of processes to be run simultaneously. The deafult is %d but this 
    should be increased by the user to the maximum number of cores available.\n""" % nThreadsDefault))
        
    print("""-h, --help
   Print this help text""")
    util.PrintCitation()       
            