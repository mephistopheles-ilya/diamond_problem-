import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot_files_to_pdf(directory, output_filename="all_plots.pdf"):
    # Get all files in the directory
    files = [f for f in os.listdir(directory) 
             if os.path.isfile(os.path.join(directory, f)) and not f.startswith('.')]
    files.sort()  # Sort files alphabetically
    
    if not files:
        print("No files found in the directory.")
        return
    
    print(f"Found {len(files)} files to plot. Creating PDF...")
    
    # Create PDF file
    with PdfPages(output_filename) as pdf:
        for filename in files:
            filepath = os.path.join(directory, filename)
            
            try:
                # Load data (handles space/comma/tab separated values)
                data = np.loadtxt(filepath)
                
                # Create figure (don't show it)
                fig = plt.figure(figsize=(8, 6))
                
                # Plot data
                if data.ndim == 1:
                    plt.plot(data, 'o-', label=filename)
                else:
                    # Assume first column is X, second is Y
                    plt.plot(data[:, 0], data[:, 1], 'o-', label=filename)
                
                plt.title(f"File: {filename}")
                plt.xlabel("X")
                plt.ylabel("Y")
                plt.legend()
                plt.grid(True)
                
                # Add this figure to the PDF
                pdf.savefig(fig)
                plt.close(fig)  # Close the figure to free memory
                
                print(f"Added plot for {filename}")
                
            except Exception as e:
                print(f"Could not plot {filename}: {str(e)}")
                continue
    
    print(f"All plots saved to {output_filename}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python plot_to_pdf.py <directory_path> [output_filename.pdf]")
        sys.exit(1)
    
    directory = sys.argv[1]
    output_filename = sys.argv[2] if len(sys.argv) > 2 else "all_plots.pdf"
    
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory")
        sys.exit(1)
    
    plot_files_to_pdf(directory, output_filename)
