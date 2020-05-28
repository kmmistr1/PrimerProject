# PrimerProject
 Primer specificity: A parameter in primer designing


This script is written to find a primer match from a multi_fasta file by trying different methods.

This script is written using sys module(sys.argv[]) to make it user friendly and to run it through terminal. It takes three arguements:

primer_fasta=sys.argv[1]  ## query seq
multi_fasta = sys.argv[2] ## test_seq
Levenshtein_distance = int(sys.argv[3])  ## number of mismatch user want to allow.

 
