import pandas as pd
import re

def shorten_taxa_names(input_file, output_file=None):
    """
    Shorten taxonomic names to their terminal (most specific) useful names.
    This is specifically for a .txt format featuretable, in this case from QIIME2
    """
    # Read the OTU table
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    
    def get_terminal_name(taxonomy_string):
        # Remove common suffixes that aren't useful
        suffixes_to_remove = [
            '_unclassified', '_uncultured', '_bacterium', '_sp.', '_strain',
            '_group', '_cluster', '_clade', '_lineage'
        ]
        
        # Split by common delimiters
        if ';' in taxonomy_string:
            parts = taxonomy_string.split(';')
            terminal = parts[-1].strip()
        else:
            terminal = taxonomy_string.strip()
        
        # Remove taxonomic prefixes (g__, s__, f__, etc.)
        terminal = re.sub(r'^[a-z]__', '', terminal)
        
        # Remove unwanted suffixes
        for suffix in suffixes_to_remove:
            terminal = terminal.replace(suffix, '')
        
        terminal = terminal.strip('_').strip()
        
        # If empty, use original
        if not terminal or terminal.isdigit() or not re.search(r'[a-zA-Z]', terminal):
            terminal = taxonomy_string
        
        return terminal
    
    # Apply the function
    df.index = df.index.map(get_terminal_name)
    
    # Handle duplicates
    index_counts = {}
    new_index = []
    
    for idx in df.index:
        if idx in index_counts:
            index_counts[idx] += 1
            new_index.append(f"{idx}_{index_counts[idx]}")
        else:
            index_counts[idx] = 0
            new_index.append(idx)
    
    df.index = new_index
    
    if output_file:
        df.to_csv(output_file, sep='\t')
        print(f"Shortened taxonomy table saved to {output_file}")
    
    return df

# Run it directly 

if __name__ == "__main__":
    # Process your file
    result = shorten_taxa_names('feature_table_full.txt', 'shortened_taxa_featuretable_L6.txt')
    print("Done!")
    
    # Show first few rows
    print("\nFirst few shortened names:")
    print(result.head())

    ##### TO DO: Tart this up to include commandline args for input and output ########