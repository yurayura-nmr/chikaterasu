import pandas as pd

# Configure pandas display settings for better output readability
pd.set_option('display.max_rows', None)  # Show all rows
pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.width', 1000)  # Prevent line wrapping in wide tables

def load_hbond_data(file):
    """
    Load hydrogen bond data from a text file and process it into a structured DataFrame.

    Parameters:
    file (str): Path to the hydrogen bond summary file.

    Returns:
    pandas.DataFrame: Processed hydrogen bond data with an additional column 'Hbond_Pair' 
                      for easy identification of unique donor-acceptor interactions.
    """
    # Read the space-separated file into a DataFrame
    df = pd.read_csv(file, delim_whitespace=True, header=None, 
                     names=["Donor", "D_Atom", "Acceptor", "A_Atom", "Exist%"])
    
    # Create a unique identifier for each hydrogen bond pair
    df["Hbond_Pair"] = df["Donor"] + "-" + df["D_Atom"] + "-" + df["Acceptor"] + "-" + df["A_Atom"]
    
    return df

# Define input file paths (update with actual paths if needed)
file1 = "summary_HBmap.dat"  
file2 = "summary_HBmap2.dat"

# Load hydrogen bond data from both simulation runs
df1 = load_hbond_data(file1)
df2 = load_hbond_data(file2)

# Merge datasets on the unique hydrogen bond pair to compare percentage existence across runs
merged_hbonds = df1.merge(df2, on="Hbond_Pair", suffixes=('_Run1', '_Run2'))

# Display the hydrogen bonds that are common between both runs with their respective percentages
print("Common Hydrogen Bonds:")
print(merged_hbonds[["Donor_Run1", "D_Atom_Run1", "Acceptor_Run1", "A_Atom_Run1", "Exist%_Run1", "Exist%_Run2"]])
