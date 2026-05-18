# Read all tsv tables in the directory and its subdirectories and if the table containt columns named  Gap_category_1 and Gap_category_2, then for each row if the Gap_category_1 is 1 and Gap_category_2 is 12, then change the values Gap_around_cut to TRUE and print Aligned_Sequence_around_cut column for this row. Finally, save the modified tables and replace the origina ones.
# %%
import os
import pandas as pd

def fix_category_1_12_in_table(file_path):
    df = pd.read_csv(file_path, sep='\t')
    
    if 'Gap_category_1' in df.columns and 'Gap_category_2' in df.columns and 'Aligned_Sequence_around_cut' in df.columns and 'Gap_around_cut' in df.columns:
        # Update Gap_around_cut to TRUE for rows where Gap_category_1 is 1 and Gap_category_2 is 12
        condition = (df['Gap_category_1'] == 1) & (df['Gap_category_2'] == 12)
        df.loc[condition, 'Gap_around_cut'] = True
        
        # Print Aligned_Sequence_around_cut for the affected rows
        print(df.loc[condition, 'Aligned_Sequence_around_cut']) 
        
        # Save the modified table back to the same file
        df.to_csv(file_path, sep='\t', index=False)

def main():
    for root, dirs, files in os.walk('/home/lefnui/CRISPRlysis_outputs/'):
        for file in files:
            if file.endswith('.tsv'):
                file_path = os.path.join(root, file)
                print(file_path)
                fix_category_1_12_in_table(file_path)

if __name__ == "__main__":
    main()
