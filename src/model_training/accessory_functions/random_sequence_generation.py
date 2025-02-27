###########################################
##Functions for generating random sequences
###########################################

#The following are functions to randomly generate sets of one hot encoded sequences that have a certain nucleotide or dinucleotide composition
#In addition, there are functions to calculate 

import random
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from sklearn.linear_model import LinearRegression
from scipy import stats

def random_seq(seqlen=1001, weights= np.array([0.25, 0.25, 0.25, 0.25])):
    
    """
    Generate a random sequence of nucleotides with specified length and weights.

    Parameters:
    seqlen (int, optional): The length of the sequence to generate. Defaults to 1001.
    weights (np.array, optional): An array of weights for the nucleotide probabilities. 
                                  The default is an equal probability for A, C, G, and T, represented as [0.25, 0.25, 0.25, 0.25].

    Returns:
    np.array: A 2D numpy array representing the generated sequence. Each row is a one-hot encoded nucleotide.
    """
    
    #one hot encodings of nucleotides, same as used in DeepSTARR's SequenceHelper.parse_alpha_to_seq

    nucleotides = [np.array([1, 0., 0, 0.]), #-->A
    np.array([0, 1, 0, 0.]), #-->C
    np.array([0, 0., 1, 0.]), #-->G
    np.array([0, 0., 0, 1])] #-->T
    
    #sample nucleotides
    seq = np.array(random.choices(nucleotides, weights = weights, k = seqlen))

    return seq

def random_X(n,seqlen=1001, weights= np.array([0.25, 0.25, 0.25, 0.25])):
    """
    Generate multiple random sequences of nucleotides with specified lengths and weights.

    Parameters:
    n (int): The number of sequences to generate.
    seqlen (int, optional): The length of each sequence. Defaults to 1001.
    weights (np.array, optional): An array of weights for the nucleotide probabilities. 
                                  The default is an equal probability for A, C, G, and T, represented as [0.25, 0.25, 0.25, 0.25].

    Returns:
    np.array: A 3D numpy array where each element is a 2D array representing a generated sequence. 
              Each row of the 2D array is a one-hot encoded nucleotide.
              Dimensionality is n x seqlen x 4
    """
    
    
    #generate n random sequences
    return np.array([random_seq(seqlen=seqlen, weights=weights) for i in range(n)])


def calculate_nucleotide_fractions(X_onehot):
    
    """
    Calculate the fractions and counts of each nucleotide (A, C, G, T) in a given one-hot encoded sequence array.

    Parameters:
    X_onehot (np.array): A 3D numpy array where each element is a 2D array representing one-hot encoded nucleotide sequences.

    Returns:
    tuple: A tuple containing two numpy arrays:
        - The first array contains the fractions of A, C, G, and T: nuc_fractions
        - The second array contains the counts of A, C, G, and T: nuc_counts
    
    nuc_fractions can serve as a weight vector for sampling with random_X.
    
    """
    
    #undo one hot encoding: looking in which position the '1' is. One hot encoded seqs are turned into sequences of 0,1,2,3.
    X_index_seqs = [np.where(X_onehot[i]==1)[1] for i in range(X_onehot.shape[0])]
    X_index_seqs = np.array([index_seq for index_seq in X_index_seqs if len(index_seq)==1001])     #exclude sequences where weird stuff happens, e.g. Ns cause abnormal lengths in this array
    
    #count overall nucleotides for normalisation
    n_nucleotides_overall = X_index_seqs.shape[0]*X_index_seqs.shape[1]
    
    #count instances of individual nucleotides
    n_As = (X_index_seqs==0).sum()
    n_Cs = (X_index_seqs==1).sum()
    n_Gs = (X_index_seqs==2).sum()
    n_Ts = (X_index_seqs==3).sum()
    
    #calculate nucleotide fractions
    A_frac = n_As/n_nucleotides_overall
    C_frac = n_Cs/n_nucleotides_overall
    G_frac = n_Gs/n_nucleotides_overall
    T_frac = n_Ts/n_nucleotides_overall
    
    nuc_fractions = np.array([A_frac,C_frac,G_frac, T_frac])
    nuc_counts = np.array([n_As,n_Cs,n_Gs,n_Ts])
    
    return nuc_fractions, nuc_counts


def calculate_dinucleotide_fractions(X_onehot):
    
    """
    Calculate the fractions and counts of each possible dinucleotide (AA, AC, AG, AT, etc.) in a given one-hot encoded sequence array.

    Parameters:
    X_onehot (np.array): A 3D numpy array where each element is a 2D array representing one-hot encoded nucleotide sequences.

    Returns:
    tuple: A tuple containing two numpy arrays:
        - The first array contains the fractions of each possible dinucleotide: dinucleotide_fracs
        - The second array contains the counts of each possible dinucleotide: dinucleotide_counts_ar
    
    dinucleotide_fracs can serve as a weight vector for sampling with random_X_dinuc.
    
    """
    
    #undo one hot encoding: looking in which position the '1' is. One hot encoded seqs are turned into sequences of 0,1,2,3.
    X_index_seqs = [np.where(X_onehot[i]==1)[1] for i in range(X_onehot.shape[0])]
    X_index_seqs = np.array([index_seq for index_seq in X_index_seqs if len(index_seq)==1001]) #exclude sequences where weird stuff happens, e.g. Ns cause abnormal lengths in this array
    

    
    # Count dinucleotide occurences, store number of occurences in dict
    # At this point, dinucleotides are encoded as tuples, e.g. (0,0) stands for AA, (1,3) stands for CG ...
    dinucleotide_counts = {}
    
    for sequence in X_index_seqs:
        # Iterate over the sequence to count dinucleotides
        for i in range(len(sequence) - 1):
            dinucleotide = (sequence[i], sequence[i + 1])
            if dinucleotide in dinucleotide_counts:
                dinucleotide_counts[dinucleotide] += 1
            else:
                dinucleotide_counts[dinucleotide] = 1
    
    #enumerate all possible dinucleotides (tuple encoding as above)
    possible_dinucs = []
    for i in range(4):
        for j in range(4):
            possible_dinucs.append((i,j))
    
    #generate vector of nucleotides of fixed length and order (as defined in the nested loop above)
    dinucleotide_counts_ar = np.array([dinucleotide_counts[dinuc] if dinuc in dinucleotide_counts.keys() else 0 for dinuc in possible_dinucs])

    dinucleotide_fracs = dinucleotide_counts_ar/dinucleotide_counts_ar.sum()
    
    return dinucleotide_fracs, dinucleotide_counts_ar


def random_seq_dinuc(seqlen=1001, weights= np.ones(16)/16):
    """
    Generate a random sequence of dinucleotides with specified length and weights.

    Parameters:
    seqlen (int, optional): The length of the sequence to generate. Defaults to 1001.
    weights (np.array, optional): An array of weights for the dinucleotide probabilities.
                                  The default is an equal probability for each of the 16 possible dinucleotides.

    Returns:
    np.array: A 2D numpy array representing the generated sequence. Each row is a one-hot encoded nucleotide.
    """
    #number of dinucleotides to sample
    n_dinucs = int(seqlen/2)
    
    nucleotides = [np.array([1, 0., 0, 0.]),
    np.array([0, 1, 0, 0.]),
    np.array([0, 0., 1, 0.]),
    np.array([0, 0., 0, 1])]

    dinucleotides = []
    for n1 in nucleotides:
        for n2 in nucleotides:
            dinuc = np.vstack([n1,n2])
            dinucleotides.append(dinuc)

    #sample dinucleotides 
    seq = random.choices(dinucleotides,weights=weights, k = n_dinucs) 
    seq = np.vstack(seq)

    #if seqlen is uneven, we need to add one single nucleotide at the end
    if seq.shape[0]!=seqlen:
        single_random_nuc = random.choice(nucleotides).reshape(1,4)
        seq = np.vstack([seq,single_random_nuc] )

    return seq

def random_X_dinuc(n,seqlen=1001, weights= np.ones(16)/16):
    """
    Generate multiple random sequences of dinucleotides with specified lengths and weights.

    Parameters:
    n (int): The number of sequences to generate.
    seqlen (int, optional): The length of each sequence. Defaults to 1001.
    weights (np.array, optional): An array of weights for the dinucleotide probabilities.
                                  The default is an equal probability for each of the 16 possible dinucleotides.

    Returns:
    np.array: A 3D numpy array where each element is a 2D array representing a generated sequence.
              Each row of the 2D array is a one-hot encoded nucleotide.
    """
    print(n)
    return np.array([random_seq_dinuc(seqlen=seqlen, weights=weights) for i in range(n)])

def testfunc():
    print("hello")

def plot_dinucfreqs(X_onehot,ax = None):
    
    """
    Plots the frequencies of all possible dinucleotides in a given one-hot encoded DNA sequence.

    Parameters:
    X_onehot (numpy.ndarray): A 3D numpy array where rows represent sequences and columns represent 
                              one-hot encoded nucleotides (shape: (num_sequences, sequence_length, 4)).
    ax (matplotlib.axes.Axes, optional): An existing Matplotlib Axes object to plot on. If None, a new 
                                         figure and axes are created. Default is None.

    Returns:
    None: The function plots the dinucleotide frequencies using Matplotlib.
    """
    
    dinuc_fracs, _ = calculate_dinucleotide_fractions(X_onehot)
    possible_dinucs = []
    
    for i in range(4):
        for j in range(4):
            possible_dinucs.append((i,j))

    #translate to letters
    nucs_dict = {0:"A", 1:"C",2:"G", 3:"T"}
    possible_dinucs_letters = [f"{nucs_dict[nuc[0]]}{nucs_dict[nuc[1]]}" for nuc in possible_dinucs]
    
    if ax!=None:
        ax.bar(possible_dinucs_letters , dinuc_fracs)
        ax.set_xlabel("Dinucleotides")
        ax.set_xlabel("Dinucleotide frequency")
        
    else:
        fig, ax = plt.subplots()
        ax.bar(possible_dinucs_letters , dinuc_fracs)
        ax.set_xlabel("Dinucleotides")
        ax.set_xlabel("Dinucleotide frequency")
        fig.show()

    
def mutate_CG(X_onehot):
    
    """
    Mutates all occurrences of the 'CG' dinucleotide in the given one-hot encoded DNA sequences to one 
    of the specified replacement dinucleotides ('TG', 'AG', 'CA', 'CT') chosen randomly.

    Parameters:
    X_onehot (numpy.ndarray): A 2D numpy array where rows represent sequences and columns represent 
                              one-hot encoded nucleotides (shape: (num_sequences, sequence_length, 4)).

    Returns:
    numpy.ndarray: A new numpy array with the same shape as X_onehot where all 'CG' dinucleotides have 
                   been mutated to one of the replacement dinucleotides.
    """
    
    target_dinuc = "CG"
    
    X_onehot_mutated = X_onehot.copy()
    
    
    nucleotides = [np.array([1, 0., 0, 0.]),
    np.array([0, 1, 0, 0.]),
    np.array([0, 0., 1, 0.]),
    np.array([0, 0., 0, 1])]

    nucleotide_names = ["A","C","G","T"]
    dinucleotide_names = []
    dinucleotide_dict = {}
    dinucleotides = []
    
    
    for n1, n1_name in zip(nucleotides, nucleotide_names):
        for n2,n2_name in zip(nucleotides,nucleotide_names):
            dinuc = np.vstack([n1,n2])
            dinucleotides.append(dinuc)
            dincu_name = f"{n1_name}{n2_name}"
            dinucleotide_names.append(dincu_name)
            dinucleotide_dict[dincu_name]=dinuc
    
    CG_array = dinucleotide_dict["CG"]
    CG_index = dinucleotide_names.index(target_dinuc)
    
    replacement_dinucleotide_names = ["TG", "AG", "CA", "CT"]
    replacement_dinucleotides = [dinucleotide_dict[dinucname] for dinucname in replacement_dinucleotide_names]
    
    #return replacement_dinucleotides
    
    for seq in X_onehot_mutated:
        #dinuc_count = 0
        for i in range(0,seq.shape[0]-1):
            #print(i)
            
            dinuc = seq[i:i+2]
            
            if (dinuc==CG_array).all():
                #dinuc_count+=1
                dinuc_replace = random.choices(replacement_dinucleotides)[0]
                seq[i:i+2]=dinuc_replace
                #print(seq[i:i+2])
        #print(dinuc_count)

                
    return X_onehot_mutated


def mutate_dinuc(X_onehot,target_dinuc = "CG", weights= np.ones(16)/16):
    ##This function does not yet work properly
    ##It can by chance create traget dinucs again
    
    
    X_onehot_mutated = X_onehot.copy()
    
    
    nucleotides = [np.array([1, 0., 0, 0.]),
    np.array([0, 1, 0, 0.]),
    np.array([0, 0., 1, 0.]),
    np.array([0, 0., 0, 1])]

    nucleotide_names = ["A","C","G","T"]
    dinucleotide_names = []
    dinucleotide_dict = {}
    dinucleotides = []
    
    
    
    for n1, n1_name in zip(nucleotides, nucleotide_names):
        for n2,n2_name in zip(nucleotides,nucleotide_names):
            dinuc = np.vstack([n1,n2])
            dinucleotides.append(dinuc)
            dincu_name = f"{n1_name}{n2_name}"
            dinucleotide_names.append(dincu_name)
            dinucleotide_dict[dincu_name]=dinuc
    
    target_dinuc_array = dinucleotide_dict[target_dinuc]
    target_index = dinucleotide_names.index(target_dinuc)
    
    
    original_target_weight =  weights[target_index]
    weights_altered = weights.copy()
    weights_altered = weights_altered + original_target_weight/15
    weights_altered[target_index] = 0
    
    
    for seq in X_onehot_mutated:
        #dinuc_count = 0
        for i in range(0,seq.shape[0]-1):
            #print(i)
            
            dinuc = seq[i:i+2]
            
            if (dinuc==target_dinuc_array).all():
                #dinuc_count+=1
                dinuc_replace = random.choices(dinucleotides,weights=weights_altered)[0]
                seq[i:i+2]=dinuc_replace
                #print(seq[i:i+2])
        #print(dinuc_count)
        break
                
    return X_onehot_mutated

def decode_onehot(X_onehot):
    """
    Decodes a onehot encoded sequence.
    Parameters: 
    X_onehot (np.array): A 2D numpy array representing a single onehot encoded DNA sequence (seqlength x 4) 
    Returns:
    decoded_seq(str): The sequence as string
    """
    
    nucleotides = {(1, 0., 0, 0.):"A",
    (0, 1, 0, 0.):"C",
    (0, 0., 1, 0.):"G",
    (0, 0., 0, 1):"T",
    (0, 0., 0, 0):"N"
    }
    
    decoded_seq=[]
    for i in range(X_onehot.shape[0]):
        nuc = nucleotides[tuple(X_onehot[i])]
        decoded_seq.append(nuc)
    decoded_seq = "".join(decoded_seq)
   
    return decoded_seq

def CG_counts_per_seq(X_onehot):
    """
    Decodes a onehot encoded sequence.
    Parameters: 
    X_onehot (np.array): A 3D numpy array representing a series of onehot encoded DNA sequences (n_seqs x seqlength x 4) 
    Returns:
    CG_counts (np.array): contains CG counts for each sequence.
    """
    
    X_letters = [decode_onehot(X_onehot[i]) for i in range(X_onehot.shape[0])]
    CG_counts = np.array([ s.count("CG") for s in X_letters])
    return CG_counts


def plot_CG_vs_accsessibility(X,y_allvar, ax=None, seqset_name = "", target_variable = "ATAC_pTE_log1p", cmap_name = "Set1"):
    """
    Create a scatterplot of accessibility against CG content for a set of DNA sequences
    Parameters:
    X (np.array): A 3D numpy array representing a series of onehot encoded DNA sequences (n_seqs x seqlength x 4) 
    y_allvar (Pandas dataframe): Data frame having at least one column named after the target variable and one column named class.
    ax: matplotlib axes object to use for plotting
    seqset_name (str): name of sequence set to be added to the title.
    """
    y=np.array(y_allvar[target_variable])
    
    cmap = plt.get_cmap(cmap_name)
    
    class_entries = y_allvar["class"].unique()
    color_dict = {cl:cmap(i) for i,cl in enumerate(class_entries)}
    class_colors = y_allvar["class"].map(color_dict)
    
    scattercolor_handles = [Line2D([0],[0],marker="o",markerfacecolor= color_dict[cl] ,color="w",label=cl) for cl in class_entries]
    
    CG_counts = CG_counts_per_seq(X)
    
    PCC = stats.pearsonr(y, CG_counts)[0]
    # Fit linear regression model
    linreg = LinearRegression()
    linreg.fit(np.array(y).reshape(-1, 1), CG_counts.reshape(-1, 1))
    slope = linreg.coef_[0]
    intercept = linreg.intercept_
    
    if ax!=None:
        ax.scatter(y, CG_counts,c=class_colors, s=0.1)
        ax.plot(y, slope * np.array(y) + intercept, color='red', label='Regression Line')
        ax.set_xlabel(target_variable)
        ax.set_ylabel("CG count per sequence")
        ax.set_title(f"{seqset_name} - Accessibility vs CG content (PCC = {PCC:.2f})")
        
        handles,labels = ax.get_legend_handles_labels()
        handles = handles+scattercolor_handles
        ax.legend(handles=handles, loc = "upper right")
   
    else:
        fig, ax = plt.subplots()
        ax.scatter(y, CG_counts,c=class_colors,s=0.1)
        ax.plot(y, slope * np.array(y) + intercept, color='red', label='Regression Line')
        ax.set_xlabel(target_variable)
        ax.set_ylabel("CG count per sequence")
        ax.set_title(f"{seqset_name} - Accessibility vs CG content (PCC = {PCC:.2f})")
        
        handles,labels = ax.get_legend_handles_labels()
        handles = handles+scattercolor_handles
        ax.legend(handles=handles, loc = "upper right")
        
        fig.show()
    

    