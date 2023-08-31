#!/usr/bin/env python3
#File: FiveEquationProgram.ipynb
#Date: 08/31/2023
#Author: Andres del Castillo
#Principal Investigator: Zachary Szpiech
#Purpose: To take equations 1-5 from the research paper "ADZE: a rarefaction approach for counting alleles private to combinations of populations", 
#   Szpiech, Jakobsson, Rosenberg (2008), and build a compute program able to compute the five equations without packaged libraries

import argparse

parser = argparse.ArgumentParser(description='Calculate ADZE Szpiech et. al (2008) Theoretical Equations')

# Files and mandatory inputs
parser.add_argument('-N', '--Nmatrix', type=str, required=True, metavar='N_Matrix_File', help='File containing the N (IxJ) matrix with rows seperated by new lines and columns by tabs')
parser.add_argument('-g', '--gsize', type=int, default=2, metavar='subsample_size', help='size of subsample [1:I-1]; default=2')
parser.add_argument('-j', '--jpop', type=int, default=2, metavar='j_population', help='j population selected [1,J]; default=2')
parser.add_argument('-k', '--kcomb', type=int, default=2, metavar='k_combination', help='k combination selected [1,J-1]; default=2')
parser.add_argument('-m', '--mcomb', type=int, default=2, metavar='m_combination', help='m-th k combination selected [1,(J choose k)]; default=2')
parser.add_argument('-i', '--alleletype', type=int, default=0, metavar='allele_type', help='Calculate Q (Equation 1) and P for a given allele type i; default=None')
args = parser.parse_args()

#Description:
#   Function to read in the integer N matrix from a text file, where row values are separated by new lines, and columns are separated by tabs
#Accepts:
#   str filename, name of the txt file containing the N matrix information  
#Returns:
#   int array matrix, N matrix to be used in further calculations

def load_matrix(filename):
    matrix = []

    with open(filename, 'r') as file:
        # Read lines one by one
        for line in file:
            # Split the line by tabs to get individual numbers
            row = line.strip().split('\t')
            # Convert strings to integers and append to the matrix
            matrix.append([int(cell) for cell in row])

    return matrix #N matrix to be used in further calculations

#Description: 
#   Function to compute the binomial coefficient between n choose k
#Accepts:
#   int n, size of set of available items to be drawn; int k, size of draw from set n
#Returns:
#   int r, value of binomial coefficient calculation

def choose(n, k):
    if k == 0: #BASE CASE 1: choosing 0 items from n
        return 1 #return 1 for the only way to choose 0 items from n
    if k > n: #BASE CASE 2: choosing more items from n than n has
        return 0 #Return 0 for not being able to choose more items than available
    r = n * choose(n - 1, k - 1) // k #Chooses first item from the set and multiplies the number of ways to choose k-1 items from the remaining n-1 items;
                                            #integer division by k to ensure we return an integer and to only count combinations by removing the arrangement factor
    return r 

#Description: 
#   Function to generate a set of possible combinations from set s with a sampling size of k
#Accepts:
#   int s, size of set with available items to be drawn; int k, size of draw from set s
#Returns:
#   array comb_arr, nested array with each element being a combination set

def generate_combinations(s, k):
    if k == 0:  #BASE CASE 1: k=0; selecting 0 items from 's'
        return [[]] #return nothing (no combinations)
    elif not s: #BASE CASE 2: s=0; selecting 'k' items from an empty list
        return [] #return nothing (no possibilities)
    else: #RECURSIVE CASE
        with_first = generate_combinations(s[1:], k-1) #Combinations w/ the first item; remove the first item from the list and find combinations with a decremented k by 1
        for combo in with_first: #Loop through found combinations w/ k-1 
            combo.append(s[0]) #Append the first element to the found combinations to make them size k
        without_first = generate_combinations(s[1:], k)#Combinations w/o the first item; remove the first item from the list and generate combinations w/ the rest of the elements
        comb_arr = with_first + without_first #concatenate results w/ and w/o first element to get final output
        return comb_arr

#Description: 
#   No Allelic Copies Probability; Equation 1, Calulates Q_{ijg}
#Params: 
#   int matrix N, an IxJ matrix with I being quantity of alleles and J the number of populations; int g, size of alleles to be subsamples; 
#       int i, allele type; int j, source population between [1,J]
#Returns:
#   float Q, probability of not finding the allele type i from the subsample g
def Q_ijg(N, g, i, j):
    Nj = sum([allele[j-1] for allele in N]) #take the sum of the column j (ensuring to map pop 1 to the first column)
    Qijg = choose(Nj - N[i-1][j-1], g) / choose(Nj, g) #calculate Q_ijg while mapping the first allele type and first population to the first entry
    return Qijg

#Description:
#   At Least One Copy Probability; Takes the complement of Q_{ijg} to find P_{ijg}
#Params: 
#   float Qijg, no copies of allele type i probability
#Returns:
#   float P, the probability of there being at least one copy of allele type i in the subsample g
def P_ijg(Qijg):
    Pijg = 1 - Qijg #Take the complete of Q_{ijg} to get P_{ijg} per the Szpiech et. al (2008) paper
    return Pijg

#Description: 
#   Estimated Allelic Richness; Equation 2, Calculates \hat{alpha}^(j)_{g}
#Params: 
#   int matrix N, an IxJ matrix with I being quantity of alleles and J the number of populations; int g, the subsample size from the source population;
#       int j, soure population between [1,J]
#Returns:
#   float alphahat_g, \hat{alpha}^(j)_{g}:the estimated allelic richess of the subsample of size g from population j
def alpha_jg(N, g, j):
    I = len(N) #number of distinct alleles in the locus 
    alphahat_g = 0 #initialize sum for \hat{alpha}^(j)_{g}
    for i in range(I): #for each allele type i in the set of distinct alleles
        Qijg = Q_ijg(N, g, i, j) #calculate Q_{ijg}
        alphahat_g += P_ijg(Qijg) #increment the sum value for \hat{alpha}^(j)_{g} accordingly 
    
    return alphahat_g

#Description:
#   Estimated Private Allelic Richness; Equation 3, calculates \hat{pi}^(j)_{g}: the estimated private allelic richness
#Accepts:
#   int matrix N, an IxJ matrix with I being quantity of alleles and J the number of populations; int g, the subsample size from the source population;
#       int j, soure population between [1,J]
#Returns:
#   float pihat_g, \hat{pi}^(j)_{g}: the estimated private allelic richness of the subsample of size g from population j   
def pi_jg(N, g, j):
    I,J = len(N), len(N[0]) #size of dinstinct allele types set, total number of populations
    pihat_gi = 0 #initialize a sum value for each allele type iteration

    JPrime_arr = list(range(1,J+1)) #Make a 1D array of the total number of populations (J)
    JPrime_arr.remove(j) #Exclude the current pop j
    for i in range(I): #for each allele type in the set of distinct alleles I
        Pijg = P_ijg(Q_ijg(N, g, i, j)) #Calulate P_{ijg} and Q_{ijg} 
        product = 1 #initialize product value for each population outside of current pop j
        for j_prime in JPrime_arr: #for each remaining population in the JPrime_arr
            Qij_primeg = Q_ijg(N, g, i, j_prime) #Calculate Q_{ij'g} for the remaining population 
            product *= Qij_primeg #increment the product value accordingly
        pihat_gi += Pijg * product #increment the sum value accordingly
    return pihat_gi

#Description: 
#   Private Allelic Combinations Across Populations; Equation 4, calculates \hat{pi}^(m)_{gk}
#Accepts:
#   int matrix N, an IxJ matrix with I being quantity of alleles and J the number of populations; int g, the subsample size from the source population;
#       int k, number of populations to be drawn from set S; int m, m-th combination of k populations
#Returns:
#   float pihat_gk, representing \hat{pi}^(m)_{gk}: the estimates number of distinct private alleles to the m-th combination of k populations drawn from subsamples of size g from each of
#       the J populations
def pi_mgk(N, g, k, m):
    I, J = len(N), len(N[0]) #size of dinstinct allele types set, number of populations

    S = list(range(1,J+1)) #initialize set S as the set of J populations labeled {1..J}
    C_k = generate_combinations(S, k) #initialize set C_k as the set of all possible combinations from S such that (size of C_k) = choose(J,k)
    C_km = C_k[m-1] #Specigic m-th combination of k populations, with '-1' to compensate for the first index of C_k correspoding to population 1
    S_without_Ckm = [j for j in S if j not in C_km] #Construct S/C_{km}, denoting the set S excluding the elements of C_{km}
    
    pihat_gk = 0 #initialize \hat{pi}^(m)_{gk} as a sum value
    for i in range(I): #for each allele type in the distinct set of alleles 
        product_Pijg = 1 #initialize the product value of P_{ijg}
        for j in C_km: #for each population in the m-th combination of k populations
            Pijg = P_ijg(Q_ijg(N, g, i, j)) #Calculate P_{ijg}
            product_Pijg *= Pijg #increment the product value of P_{ijg} accordingly

        product_Qijg = 1 #initialize the product value of Q_{ij'g}
        for j_prime in S_without_Ckm: #for each remaining population in S not in C_km
            Qij_primeg = Q_ijg(N, g, i, j_prime) #calculate Q_{ij'g}
            product_Qijg *= Qij_primeg #increment the product value of Q_{ij'g} accordingly

        pihat_gk += product_Pijg * product_Qijg #Take the product of the cumulative products of P_{ijg} and Q_{ij'g}, and increment the \hat{pi}^(m)_{gk} sum value accordingly
    
    return pihat_gk

#Description:
#   Missing Allelic Richness; Equation 5, gives the measure of the number of distinct alleles found in all populations but j
#Accepts:
#int matrix N, an IxJ matrix with I being quantity of alleles and J the number of populations; int g, the subsample size from the source population;
#   int j, source population between [1,J]
#Returns:
#   float muhat_g, representing \hat{mu}^(j)_{g}: the sample size-corrected measure of the number of distinct alles found in all pops but population j      
def mu_jg(N, g, j):
    I, J = len(N), len(N[0]) #size of dinstinct allele types set, number of populations

    JPrime_arr = list(range(1,J+1)) #Make a 1D array of the total number of populations (J)
    JPrime_arr.remove(j) #Exclude the current pop j
    
    muhat_g = 0 #initialize the sum value for \hat{mu}^(j)_{g}
    for i in range(I): #for each distinct allele in the locus
        Qijg = Q_ijg(N, g, i, j) #Calculate Q_{ijg}
        product = 1 #Initialize the product for P_{ij'g}
        for j_prime in JPrime_arr: #for each remaining population
            Pij_primeg = P_ijg(Q_ijg(N, g, i, j_prime)) #Calculate P_{ij'g}
            product *= Pij_primeg #Incremenet the product for P_{ij'g} accordingly
        muhat_g += Qijg * product #increment the sum value for \hat{mu}^(j)_{g} by multiplying q_{ijg} by the iterated product of P_{ij'g}
    return muhat_g

def run_program(N,g,k,j,m):
    return [alpha_jg(N, g, j), pi_jg(N, g, j), pi_mgk(N, g, k, m), mu_jg(N, g, j)]

if __name__ == '__main__':
    N = load_matrix(args.Nmatrix)
    outputted_equations = run_program(N, args.gsize, args.kcomb, args.jpop, args.mcomb)
    
    if args.alleletype>0:
        Q = Q_ijg(N, args.gsize, args.alleletype, args.jpop)
        print(f'Equation 1; No Copies Probability:{Q}')
        
    print(f'Equation 2; Estimated Allelic Richness of g: {outputted_equations[0]}')
    print(f'Equation 3; Estimated Private Allelic Richness of g: {outputted_equations[1]}')
    print(f"Equation 4; Estimated # of distinct private alleles to m-th combination of 'k' populations of size 'g': {outputted_equations[2]}")
    print(f"Equation 5; Size-corrected measure of # of distinct alleles found in all pops but pop j: {outputted_equations[3]}")
