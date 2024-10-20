import re
import os

class AntigenExtractor:
    def __init__(self):
        self.file_contents = {}
        self.antigen_pattern = r'\b(?:CD\d+|HLA-\w+|[A-Z]{2,}\d*)\b'
    
    def read_files_from_directory(self, directory_path):
        """
        Read all text files from a given directory and store their contents.
        
        Parameters:
        directory_path (str): Path to the directory containing the text files.
        """
        for filename in os.listdir(directory_path):
            if filename.endswith(".txt"):
                file_path = os.path.join(directory_path, filename)
                with open(file_path, 'r', encoding='utf-8') as file:
                    self.file_contents[filename] = file.read()
    
    def extract_antigens(self):
        """
        Extract all unique antigen names mentioned in the files.
        
        Returns:
        dict: A dictionary where keys are file names and values are sets of unique antigen names.
        """
        antigen_data = {}
        for file_name, content in self.file_contents.items():
            antigens = set(re.findall(self.antigen_pattern, content))
            antigen_data[file_name] = antigens
        return antigen_data
    
    def query_unique_antigens(self):
        """
        Get the unique antigen names across all files.
        
        Returns:
        set: A set of unique antigen names found across all files.
        """
        all_antigens = set()
        antigen_data = self.extract_antigens()
        for antigens in antigen_data.values():
            all_antigens.update(antigens)
        return all_antigens


# Usage example

# Create an instance of the extractor
extractor = AntigenExtractor()

# Provide the path to the directory where your text files are stored
directory_path = 'downloaded_articles/txt'

# Read the files from the directory
extractor.read_files_from_directory(directory_path)

# Extract unique antigens from each file
antigen_data = extractor.extract_antigens()
print("Antigens found in each file:")
for file, antigens in antigen_data.items():
    print(f"{file}: {antigens}")

# Query unique antigens across all files
unique_antigens = extractor.query_unique_antigens()
print("\nUnique antigens across all files:")
print(unique_antigens)