import os
import sys
import logging
from PyPDF2 import PdfReader

def setup_logging():
    """
    Sets up logging configuration.
    Logs are saved to 'pdf_to_text.log' in the script's directory.
    """
    logging.basicConfig(
        filename='pdf_to_text.log',
        filemode='a',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO
    )

def convert_pdf_to_text(pdf_path, txt_path):
    """
    Converts a single PDF file to a text file.

    Args:
        pdf_path (str): Path to the PDF file.
        txt_path (str): Path where the text file will be saved.
    """
    try:
        reader = PdfReader(pdf_path)
        text = ""
        for page_num, page in enumerate(reader.pages, start=1):
            page_text = page.extract_text()
            if page_text:
                text += page_text
            else:
                logging.warning(f"No text found on page {page_num} of {pdf_path}")
        
        with open(txt_path, 'w', encoding='utf-8') as f:
            f.write(text)
        
        logging.info(f"Successfully converted '{pdf_path}' to '{txt_path}'")
    except Exception as e:
        logging.error(f"Error converting '{pdf_path}': {e}")

def process_directory(directory):
    """
    Processes all PDF files in the given directory.

    Args:
        directory (str): Path to the directory to process.
    """
    if not os.path.isdir(directory):
        logging.error(f"The directory '{directory}' does not exist or is not a directory.")
        print(f"Error: The directory '{directory}' does not exist or is not a directory.")
        return

    # List all files in the directory
    files = os.listdir(directory)
    pdf_files = [file for file in files if file.lower().endswith('.pdf')]

    if not pdf_files:
        logging.info(f"No PDF files found in the directory '{directory}'.")
        print(f"No PDF files found in the directory '{directory}'.")
        return

    print(f"Found {len(pdf_files)} PDF file(s) in '{directory}'. Processing...")

    for pdf_file in pdf_files:
        pdf_path = os.path.join(directory, pdf_file)
        base_name = os.path.splitext(pdf_file)[0]
        txt_file = base_name + '.txt'
        txt_path = os.path.join(directory, txt_file)
        
        convert_pdf_to_text(pdf_path, txt_path)
    
    print("Processing completed. Check 'pdf_to_text.log' for details.")

def main():
    setup_logging()

    # Check if directory path is provided as a command-line argument
    if len(sys.argv) > 1:
        directory = sys.argv[1]
    else:
        # Prompt the user to enter the directory path
        directory = input("Enter the path to the directory containing PDF files: ").strip()
    
    process_directory(directory)

if __name__ == "__main__":
    main()