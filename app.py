import tkinter as tk
from tkinter import ttk, messagebox
import sys
import os
from io import StringIO
import main

class SequenceAlignmentApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Sequence Alignment Tool")
        self.root.geometry("700x550")
        self.root.configure(bg="#f0f0f0")
        
        # Style
        style = ttk.Style()
        style.configure("TLabel", background="#f0f0f0", font=("Helvetica", 10))
        style.configure("TButton", font=("Helvetica", 10, "bold"))
        
        # File Input Frame
        frame_input = ttk.Frame(root, padding=10)
        frame_input.pack(fill="x")
        ttk.Label(frame_input, text="Enter FASTA files (comma-separated):").pack(side="left")
        self.entry_files = ttk.Entry(frame_input, width=40)
        self.entry_files.pack(side="left", padx=5)
        self.entry_files.insert(0, "sequence.fasta")
        
        # Checkboxes for Options Frame
        frame_options = ttk.Frame(root, padding=10)
        frame_options.pack(fill="x")
        ttk.Label(frame_options, text="Select Analysis Methods:").pack(anchor="w")
        
        # Options correspond to main.py choices
        self.options = [
            "Needleman-Wunsch (Global Alignment)",
            "Smith-Waterman (Local Alignment)",
            "Star Alignment (Multiple Sequence Alignment)",
            "PSSM Analysis",
            "Semi-global Alignment",
            "Bayesian Conservation Analysis",
            "HMM Training"
        ]
        
        self.checkboxes = []
        for i, opt in enumerate(self.options):
            var = tk.BooleanVar()
            chk = ttk.Checkbutton(frame_options, text=opt, variable=var)
            chk.pack(anchor="w", padx=5)
            self.checkboxes.append(var)
        
        # Add "All" checkbox
        self.all_var = tk.BooleanVar()
        chk_all = ttk.Checkbutton(frame_options, text="All of the above", 
                                  variable=self.all_var, 
                                  command=self.toggle_all)
        chk_all.pack(anchor="w", padx=5)
        
        # Run Button
        ttk.Button(root, text="Run Analysis", command=self.run_analysis).pack(pady=10)
        
        # Output Text Box Frame
        frame_output = ttk.Frame(root, padding=10)
        frame_output.pack(fill="both", expand=True)
        self.text_output = tk.Text(frame_output, height=20, width=80, 
                                  font=("Courier", 10), bg="#ffffff", 
                                  relief="flat", borderwidth=1)
        self.text_output.pack(fill="both", expand=True)
        
        # Add scrollbar to text output
        scrollbar = ttk.Scrollbar(frame_output, orient="vertical", 
                                 command=self.text_output.yview)
        scrollbar.pack(side="right", fill="y")
        self.text_output.config(yscrollcommand=scrollbar.set)
    
    def toggle_all(self):
        """Toggle all checkboxes based on "All" checkbox state"""
        state = self.all_var.get()
        for var in self.checkboxes:
            var.set(state)
    
    def run_analysis(self):
        """Call main.py's functionality with GUI inputs"""
        # Clear output first
        self.text_output.delete(1.0, tk.END)
        
        # Get FASTA files from entry
        fasta_input = self.entry_files.get().strip()
        if not fasta_input:
            self.text_output.insert(tk.END, "Error: Please enter at least one FASTA file.")
            return
        
        # Get selected analysis methods
        choices = []
        for i, var in enumerate(self.checkboxes):
            if var.get():
                choices.append(i+1)  # main.py uses 1-based indexing for choices
        
        # Also check if "All" is selected
        if self.all_var.get():
            choices.append(8)  # 8 is the "All" option in main.py
        
        if not choices:
            self.text_output.insert(tk.END, "Please select at least one analysis method.")
            return
        
        # Capture console output from main.py
        original_stdout = sys.stdout
        captured_output = StringIO()
        sys.stdout = captured_output
        
        try:
            # Save original input function and mock it
            original_input = __builtins__.input
            
            # Mock input function to provide our GUI values to main.py
            def mock_input(prompt):
                if "Enter sequence file names" in prompt:
                    return fasta_input
                elif "Enter your choice" in prompt:
                    return ",".join(map(str, choices))
                return ""
            
            __builtins__.input = mock_input
            
            # Call main function
            try:
                main.main()
            except SystemExit:
                pass  # Ignore SystemExit from main.py
            except Exception as e:
                self.text_output.insert(tk.END, f"Error executing analysis: {str(e)}\n")
                import traceback
                self.text_output.insert(tk.END, traceback.format_exc())
            
            # Get captured output and display in text box
            output = captured_output.getvalue()
            self.text_output.insert(tk.END, output)
            
        finally:
            # Restore stdout and input
            sys.stdout = original_stdout
            __builtins__.input = original_input
            
            # Ensure we always scroll to the start
            self.text_output.see("1.0")

if __name__ == "__main__":
    root = tk.Tk()
    app = SequenceAlignmentApp(root)
    root.mainloop()