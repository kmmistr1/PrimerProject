##   This script contains codes to solve the Bioinformatic problem (primer designing: finding a primer match)

## This script inculdes the codes written to define functions to perform specific task
## and these functions are applied/called on the files to tackle the problem

## This script intakes 3 arguements in the form of sys.argv[] to make user friendly where usercan use their file of choice and number of mismatches user want to look for.
###These are explained in detail below.

## I. Step to import modules Necessary modules that we will need
### To import these modules:
import sys                                                    ### sys module allows us to run the code through command line
import re                                                     ### to import the module regular expression
import fuzzysearch                                            ### external modue. 1st installed in terminal and then imported here
from fuzzysearch import find_near_matches                     ### to import find_near_matches function from fuzzy search module.
from Bio import pairwise2                                     ### Importing pairwise2 module from biopython to perform alignment
from Bio.pairwise2 import format_alignment                    ### from pairwise2 import the function format alignment




## II. FILE HANDLING using "sys" module and loop method
### This script uses three arguments (sys.argv[]) to make the script user friendly
### Before running the script, user need to input these 3 arguments
Primer = sys.argv[1]                                        ### the primer (query) file that need to be imported using user pathway
multi_fasta = sys.argv[2]                                   ###  the multi_fasta (test) file that need to be imported using user pathway
Levenshtein_dist = int(sys.argv[3])                         ###  this is the number of mismatches user want to allow

## A). QUERY =  Primer fasta file
with open(Primer,"r")as file:                               ### open and read Primer file using "with open()"
    for lines in file:                                      ###  iterates over lines in file
        if not lines.startswith(">"):                       ### considering the lines that does not start with ">"
            primer=lines.upper().strip()                    ### calling the lines "primer" and converting them in to upper case


### B). Test = multi fasta file
### converting it into a dictionary
test_dict={}                                                ### creating an empty dictionary
with open(multi_fasta, "r")as file:                         ### open and read test_sequence file using "with open()"
    for lines in file:                                      ### iterates over lines in file
        if lines.startswith(">"):                           ### considering the lines that starts with ">" (condition applied)
            header = lines.strip()                          ### calling the lines header

        else:                                               ### if not above then do this
            seq = lines.upper().strip()                     ### calling the lines seq
            test_dict[header] = str(seq)                    ### adding keys=header and values=seq to an empty dictionary




## III. FUNCTIONS to convert primer seq into different forms
### Match could be in the form of complement bases too, so, to define them:
###  A. function to complementary bases for sequences.
def complement(sequence):                                  ### defining a function
    '''returning complement bases of the given sequence'''
    List=[]                                                ### creating an empty list
    seq_dict={"A":"T","T":"A","G":"C","C":"G"}             ### dictionary of complement bases
    for x in sequence:                                     ### for bases in sequence, iterates over each base
        if x in seq_dict:                                  ### for complement bases in dictionary
            List.append(seq_dict[x])                       ### if bases of sequence are present in dictionary then add them to the list
    return ''.join(List)                                   ### returning the list

### B. defining a function to covert complement into its reverse seq
def revcom(sequence):                                      ### difining a function
    '''this function converts the given sequence into their
    complement bases and then return them in reverse order '''
    return complement(sequence[::-1])                      ### return the complement bases into reverse order
### So now we have 3 query sequences to find match for: pimer, complement(primer) and revcom(primer).



## IV. METHODS apply to find the match/matches:
### Below are the different techniques that were tried to find all the possible matches

### A. Basic method: using "if...else"
def find_string(Dict, string):                             ### function is defined to find a pattern in the string
    '''return the list of the lines if the string exist in the lines of dictionary,
    if not then will print no match found'''
    match_list=[]                                          ### creating an empty list
    for k,v in Dict.items():                               ### for keys and values in the dictionary
        if string in v:                                    ### if the pattern is == to values/in the lines in the dictionary
            match_list.append([k,v])                       ### then add those to the empty list along with theie headers
        else:                                              ### if not pattern in the lines
            match_list.append([k,"No match found"])        ### then print this
    return(match_list)                                     ### returning the list

## Applying "find_string()" function on all 3 query sequences:
print("Three diffrent methods tried to find the matches are:","\n")
print("A. Basic method ('if...else'):","\n")
print("Primer sequence:")
print (find_string(test_dict, primer),"\n")                ### function applied to primer seq
print("Complement(primer) sequence:")
print (find_string(test_dict, complement(primer)),"\n")    ### function applied to complement(primer) seq
print("Revcom(primer) sequence:")
print (find_string(test_dict, revcom(primer)),"\n")        ### function applied to revcom(primer) seq

### Unfortunately, this approach won't give the match string or it's location
### but it's good way to get an idea in which part of the test file our match can exist


### B. RegEx (regular expression):
### Here, is the more sophisticated method to find the match as well as it's location if the match exist in the string
def find_match_regex(Dict, string):                         ### defining a function to find a match using re
    '''return the listof the match objects if the match is found,
     if there is no match it will print "no match found"'''
    match_object=[]                                        ### create an empty list
    for k,v in Dict.items():                               ###for keys(headers) and values(seqs) in the dictionary
        R=re.search(string,v)                              ### using search function from "re" to find a pattern in the string
        if R==None:                                        ### if there is no match
            match_object.append([k,"No match found"])      ### add this statement to the dictionary
        else:                                              ### if match found
            match_object.append([k,R])                     ### add them to dictionary too along with the header
    return(match_object)                                   ### return the list

## Applying "find_match_regex()" on all the 3 queries and test file
print("B. Regex method:","\n")
print("Primer: ")
print(find_match_regex(test_dict, primer),"\n")             ### applied to primer seq
print("Complement(primer): ")
print(find_match_regex(test_dict, complement(primer)),"\n") ### applied to complement(primer)
print("Revcom(primer): ")
print(find_match_regex(test_dict, revcom(primer)),"\n")    ### applied to revcom(primer)

### This method will give the 100% match.


### C. Fuzzy searching:
### it's time to get matches with 10% mistaches
### I have used fuzzy search module to do fuzzy searching which allows the mismatches to find the match
def fuzzy_match(pattern, Dict, dist_lev):                  ### defining a function to find a fuzzy match
    '''it finds the pattern match from the dictionary by allowing some mismatches,
     here we will use Levenshtein distance as a parameter'''
    List=[]                                                ### create empty dict
    for k,v in Dict.items():                               ### for items in dictionary
        L=(find_near_matches(pattern,v, max_l_dist = dist_lev))     ### using find_near_matches() function which use levenshtein distance
                                                                    ### here we can give the number of mismatches that we want to allow
        if L:                                             ### if matches are found
            List.append([k,L])                            ### append them to the list along with keys
        else: List.append([k,"No match"])                 ### if not then add "No match"
    return(List)                                          ### returning the list

## calling above function on all the 3 queries ###
## Levenshtein_distance is the argument that user can put of his/her choice
## creating empty lists for respective query sequences

comp_list=[]                                              ### for primer
revcom_list=[]                                            ### for complement(primer)
normal_list=[]                                            ### for revcom(primer)
print("C. Fuzzy searching approach: ","\n")

print("Primer:")                                          ### print statement
for matches in fuzzy_match(primer, test_dict, Levenshtein_dist):   ### calling above function on primer seq
    for y in matches:                                     ### itertaes over each line
        normal_list.append(y)                             ### appending lines to the list
print(normal_list,"\n")                                   ### printing list
print("Complement(primer):")
for matches in fuzzy_match(complement(primer), test_dict, Levenshtein_dist):  ### calling function on complement(primer)
    for y in matches:                                    ### iterates over line
        comp_list.append(y)                              ### appending lines to the list
print(comp_list,"\n")                                    ### printing the list
print("\n","Revcom(primer):")
for matches in fuzzy_match(revcom(primer), test_dict, Levenshtein_dist):  ###calling above function on revcom(primer) seq
    for y in matches:                                    ### iterates over each line
        revcom_list.append(y)                            ### appending the lines to the list
print(revcom_list,"\n")
print("\n","Fuzzy search outputs are chosen over basic method and RegEx for futher processing")                             ### printing the list

### I chose this method over the basic and RegEx method for further investigation
### As this approach gives the ouptput fulfilling all conditions in a easy way




## V. ALIGNMENT:
## to check the similarity between the matched query sequences and the matched strings
## I decided to  pairewise2 alignment (localms) by penalizing  matches, mismatches and open/extend gaps

### A. function is defined to extract the match objects and use them for alignments
def alignment_pairwise2(fuzzy_list,string, Dict):                                           ### defining function for alignment
    '''perform the alignment of the string and pattern
    by using extracted match objects from the list'''
    empty_list=[]                                                                           ### create empty list
    for x in range(1, len(fuzzy_list),2):                                                   ### iterates over the matched list in the given range
        start_ = int(str(fuzzy_list[x]).split(",")[0].split("=")[1])                        ### extracting the start position using "split" method
        end_ = int(str(fuzzy_list[x]).split(",")[1].split("=")[1])                          ### to extract the end position using "split"
        match_ = str(fuzzy_list[x]).split(",")[3].split("=")[1].strip("'").split("'")[0]    ### extract matched string

        for k,v in Dict.items():                                                            ### for keys and values in test_dict
            if match_ in v:                                                                 ### if matched string is present in values
                for a in pairwise2.align.localms(v[start_:end_:1],string,1, 0, -1, -0.5, one_alignment_only=True):
                    ### performing alignment and iterates over each alignment
                    empty_list.append(k +"\n"+format_alignment(*a))                         ### appending the outputs to the list in a precise format

    return(empty_list)                                                                      ### returning the list

### B. defining a function in conjunction with "alignment_pairwise2()" which will be applied on query seqs
def align_output(query_list, query_string,dictionary):
    '''return alignments, if no matches are found then will print "No alignments"'''
    if not "No match" in query_list:                                                        ### If the list dont have this statement
        for aligns in alignment_pairwise2(query_list,query_string,dictionary):              ### then peform alignment using predefined function
            print(aligns)                                                                   ### print alignments
    else:                                                                                   ### if have "No match"
        for items in query_list:                                                            ### iterates over each line
            if ">" in items:                                                                ### if line have ">"
                print("{0}:\n No alignments\n".format(items))                               ### then print this

## C. printing alignments for the matches using "out_put()""
print("The alignment below shows the similarity and error type between the query and the matched string along with the alignment score.","\n"
print("Pairwise2 alignments :","\n" )
print("Pimer:")
align_output(normal_list, primer, test_dict)                                                ### alignment for primer seq
print("\n")
print("Complement of Pimer:")
align_output(comp_list, complement(primer), test_dict)                                      ### alignment for complement(primer) seq
print("\n")
print("Revcom(Pimer):")
align_output(revcom_list, revcom(primer),test_dict)                                         ### alignment for revcom(primer) seq
