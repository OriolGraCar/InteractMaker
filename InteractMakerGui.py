import tkinter as tk
from PIL import ImageTk, Image
from tkinter import filedialog
import re
import sys
import os
import copy
import threading
import queue
import time
import subprocess


class Writer(object):
    """This Class is used to override and redirect the stdout and stderr stream to the Tk app"""
    def __init__(self, parent, channel):
        self.parent = parent
        self.channel = channel

    def write(self, data):
        # Function that is called in the scripts when writing to the channels
        self.parent.console.update_text(self.channel, data)


class MenuBar(tk.Menu):
        """Widget of the Menu bar of the GUI"""
        def __init__(self, parent):
            tk.Menu.__init__(self, parent)
            self.parent = parent
            # Define Cascade SubMenus
            self.filemenu = tk.Menu(self, tearoff=0)
            self.modulemenu = tk.Menu(self, tearoff=0)
            # Options for Filemenu
            self.filemenu.add_command(label="Open", command=self.openfiles)
            self.filemenu.add_command(label="Open Folder", command=self.opendir)
            self.filemenu.add_command(label="Save", command=self.save)
            self.filemenu.add_separator()
            self.filemenu.add_command(label="Exit", command=self.parent.quit)
            # Options for Modulemenu
            self.modulemenu.add_command(label="Interact Maker", command=lambda x="Interact Maker": self.change_module(x))
            self.modulemenu.add_command(label="PDB splitter", command=lambda x="PDB splitter": self.change_module(x))
            self.modulemenu.add_command(label="PDB validator", command=lambda x="PDB validator": self.change_module(x))
            # Display
            self.add_cascade(label="File", menu=self.filemenu)
            self.add_cascade(label="Modules", menu=self.modulemenu)
            self.add_command(label="Help")

        def save(self):
            """Function to save output into file (provisional)"""
            out = filedialog.asksaveasfile(title="Save Structure as:", mode='w', defaultextension=".pdb")
            if out is not None:
                if self.parent.current_module_name == "Interact Maker":
                    filename = out.name
                    if self.parent.current_module.output is not None:
                        self.parent.current_module.output.save_to_file(filename)
                out.close()



        def change_module(self,modul):
            """
            Function that changes the modul frame that is beign displayed.
            Changing module sis fats because we have everything loaded into memory
            """
            father_modules = {"Interact Maker": self.parent.interact_frame,
                              "PDB splitter": self.parent.pdbsplit_frame,
                              "PDB validator": self.parent.pdbvalidate_frame}
            if self.parent.current_module_name != modul:
                """
                because the method pack_forget() removes all the configuration,
                we have to unpack everything and repack it
                Packing and unpacking is fast because we have all the frames and their children widgets in memory
                """
                father_modules[self.parent.current_module_name].pack_forget()
                self.parent.output_options_frame.pack_forget()
                self.parent.output_frame.pack_forget()
                father_modules[modul].pack(fill=tk.X)
                self.parent.output_options_frame.pack(fill=tk.X)
                self.parent.output_frame.pack(fill='both', expand=True)
                self.parent.current_module_name = modul

        def openfiles(self):
            """ Funnction that open pdb files, it have different behaviour with each modul"""
            filename = filedialog.askopenfilename(title="Select PDB file", filetypes=(("Pdb files", "*.pdb"),))
            if filename:
                name = re.search('^.*/(.*)\.pdb', filename)
                if name:
                    if self.parent.current_module_name == "Interact Maker":
                        protein = PS(name.group(1), filename)
                        # Here we call the loading funciton of the correspondent module
                        self.parent.current_module.update_proteins([(protein, filename)])
                else:
                    # This condition should never happen since the files are masked
                    print("Unexpected ERROR when opening file, '^.*/(.*)\.pdb' didn't matched")

        def opendir(self):
            """Function that allows us to open a whole directory, loading all the pdbs into memory"""
            folder_name = filedialog.askdirectory(title="Select PDB folder")
            proteins_to_send = []
            if folder_name:
                for file in os.listdir(folder_name):
                    name = re.search('(.*)\.pdb', file)
                    if name:
                        path = folder_name + "/" + name.group(1) + ".pdb"
                        protein = PS(name.group(1), path)
                        proteins_to_send.append((protein, path))
                # Here we call the function to update the files of the correspondent module
                self.parent.current_module.update_proteins(proteins_to_send)


class MatchModule(tk.Frame):
    """Frame of the module that reconstructs the protein"""
    def __init__(self, parent, master):
        # Basic setting
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.master = master
        # Variables needed for threading
        # The program sends child threads that the user can stop, relaunch and kill. For this reason
        # we need several booleans to check the current state of the job, and to avoid launching multiple threads.
        self.restart = True
        self.running = False
        self.job = None
        # Widgets declaration
        self.load_butt = tk.Button(self, text='Load', command=self.master.menu.openfiles)
        self.run_butt = tk.Button(self, text='Run ', command=self.run_v2)
        self.stop_butt = tk.Button(self, text='Kill', command=self.kill_job_v2)
        self.erase_butt = tk.Button(self, text='Erase', command=self.remove_protein)
        self.info_butt = tk.Button(self, text='Info', command=self.someinfo)
        self.image_butt = tk.Button(self, text='Image', command=self.update_img)
        self.pymol_butt = tk.Button(self, text='Open in Pymol')
        self.info_label = tk.Label(self, text='Select a option to apply to: None', width=45, anchor=tk.W)
        self.erase_all_butt = tk.Button(self, text="Erase All", command=self.remove_all)
        self.log_label = tk.Label(self)
        self.protein_display = tk.Listbox(self, height=27)
        # Variables that will contain the information that the user loads and generates
        self.proteins = {}
        self.protein_list = []
        self.output = None
        # Widgets placement after declaration
        self.load_butt.grid(row=0, column=0)
        self.run_butt.grid(row=0, column=1)
        self.stop_butt.grid(row=0, column=2)
        self.protein_display.grid(row=1, column=0, columnspan=3, rowspan=20)
        self.info_label.grid(row=1, column=4, sticky=tk.W)
        self.erase_butt.grid(row=3, column=4, sticky=tk.W)
        self.info_butt.grid(row=4, column=4, sticky=tk.W)
        self.image_butt.grid(row=5, column=4, sticky=tk.W)
        self.pymol_butt.grid(row=6, column=4, sticky=tk.W)
        self.erase_all_butt.grid(row=7, column=4, sticky=tk.W)
        self.log_label.grid(row=8, column=4, sticky=tk.W)
        # Some functionality
        self.protein_display.bind('<<ListboxSelect>>', self.update_info)

    def update_img(self):
        """Function that changes the current image"""
        prot_indx = self.protein_display.curselection()
        if prot_indx:
            prot_indx = int(prot_indx[0])
            file_name = self.proteins[self.protein_list[prot_indx]]["path"]
            self.master.img_display.start_pymol_image_thread(file_name)

    def update_info(self, *args):
        """Function that updates the info label when a item of the listbox is selected"""
        # the reason why i use *args is because the method bind sends as argument the trigger event so if it
        prot_indx = self.protein_display.curselection()
        if prot_indx:
            prot_indx = int(prot_indx[0])
            self.info_label.config(text='Select a option to apply to: %s' % self.protein_list[prot_indx])
        else:
            self.info_label.config(text='Select a option to apply to: None')

    def update_proteins(self, proteins_to_add):
        """Function that loads a list of proteins to the memory."""
        for protein in proteins_to_add:
            if protein[0].get_id() not in self.proteins:
                self.protein_display.insert(tk.END, protein[0].get_id())
                self.proteins[protein[0].get_id()] = {"prot":protein[0], "path": protein[1]}
                self.protein_list.append(protein[0].get_id())
            if protein[0].get_id() == "MacroComplex":
                self.proteins[protein[0].get_id()] = {"prot":protein[0], "path": protein[1]}

    def remove_protein(self):
        """Function that removes the selected protein"""
        prot_indx = self.protein_display.curselection()
        if prot_indx:
            prot_indx = int(prot_indx[0])
            self.protein_display.delete(prot_indx)
            self.proteins.pop(self.protein_list[prot_indx])
            sys.stdout.write("Protein %s has been removed\n" % self.protein_list[prot_indx])
            self.log_label['text'] = "Protein %s has been removed\n" % self.protein_list[prot_indx]
            self.protein_list.pop(prot_indx)

    def remove_all(self):
        """Function that removes all the loaded information"""
        self.protein_display.delete(0, tk.END)
        sys.stdout.write("All proteins had been removed\n")
        self.log_label['text'] = "All proteins had been removed\n"
        self.proteins = {}
        self.protein_list = []

    def someinfo(self):
        """Function that displays into the console the basic info of the selected listbox entry"""
        prot_indx = self.protein_display.curselection()
        if prot_indx:
            prot_indx = int(prot_indx[0])
            prot_obj = self.proteins[self.protein_list[prot_indx]]["prot"]
            prot_name = prot_obj.get_id()
            self.log_label['text'] = "%s Info printed in the info tap" % prot_name
            info_console.write("The Protein %s has:\n%s Chains\n" % (prot_name, len(prot_obj)))
            total_residues = 0
            total_atoms = 0
            for chain in prot_obj:
                atom_counter = 0
                total_residues += len(chain)
                info_console.write("\tChain %s has:\n\t\t%s Residues" % (chain.get_id(), len(chain)))
                for residue in chain:
                    atom_counter += len(residue)
                    total_atoms += len(residue)
                info_console.write(" and %s Atoms\n" % atom_counter)
            info_console.write("Total:\n%s Residues and %s Atoms\n" % (total_residues, total_atoms))

    def run_v2(self):
        """
         Function that launches a new child thread if no threads are running,
         stops the current thread and relaunches stopped threads.
         it also checks that all the info loade dinto memory is correct
        """
        all_correct = True
        if not self.running:
            if self.restart:
                for protein in self.proteins:
                    if len(self.proteins[protein]["prot"]) != 2:
                        all_correct = False
                        self.log_label['text'] = "All proteins should have two chains"
                        sys.stderr.write("All proteins should have two chains\n")
            if len(self.protein_list) == 0:
                all_correct = False
                self.log_label['text'] = "No proteins Found"
                sys.stderr.write("No proteins Found\n")
            if all_correct:
                self.run_butt['text'] = "Stop"
                self.log_label['text'] = "Job Running"
                sys.stdout.write("Job Running\n")
                self.running = True
                if self.restart:
                    self.job = ReconstructComplex(self)
                    self.restart = False
                    self.job.pdb_list = [protein["prot"] for pid, protein in self.proteins.items()]
                    self.job.homo_chains = cmake.good_chain_names(self.job.pdb_list)
                    self.job.interactions_dict = cmake.all_interactions(self.job.pdb_list, self.job.homo_chains)
                    self.job.new_pdb = copy.deepcopy(self.job.pdb_list[0])
                    self.job.new_pdb.id = "MacroComplex"
                    self.job.chain_id_dict = {self.job.new_pdb.childs[0].get_id(): self.job.new_pdb.childs[0].get_id(),
                                              self.job.new_pdb.childs[1].get_id(): self.job.new_pdb.childs[1].get_id()}
                    self.job.start()
        else:
            self.run_butt['text'] = "Stop"
            self.log_label['text'] = "Job Stopped"
            self.running = False
            sys.stdout.write("Job Stopped\n")
            self.run_butt['text'] = "Run "

    def kill_job_v2(self):
        """Restart everything back to normal."""
        if not self.restart:
            self.restart = True
            self.running = False
            self.run_butt['text'] = "Run "
            self.remove_all()
            sys.stdout.write("Job Killed\n")


class ReconstructComplex(threading.Thread):
    """
    Class that handles the Recontruction job.
    The reason to use threads is to avoid freezing the mainloop of the gui.
    """
    def __init__(self, master):
        threading.Thread.__init__(self)
        self.master = master
        self.setDaemon(True)
        # Variables needed for the reconstruct
        self.pdb_list = list()
        self.homo_chains = {}
        self.new_pdb = None
        self.interactions_dict = {}
        self.chain_id_dict = {}
        self.tmp_count = 0
        self.completed_borders = 0
        self.completed_chains = 0
        self.confirmed_residues = 0
        self.lax_residues = 0
        self.current_chain_id = None
        self.current_border_list = []
        self.interactions_checked = False


    def run(self):
        """
        Recursive function that rebuilds the macrocmplex.
        it checks two labels: self.master.running that contains information if the job should be running or be waiting
        and self.master.restart: that allows us to finish the recursive calls
        """
        if not self.master.restart:
            if self.master.running:
                chain = self.new_pdb.childs[self.completed_chains]
                if chain.get_id() != self.current_chain_id:
                    self.current_chain_id = chain.get_id()
                    self.current_border_list = list(self.interactions_dict[self.chain_id_dict[chain.get_id()]].keys())
                border = self.current_border_list[self.completed_borders]
                if self.interactions_checked:
                    if self.confirmed_residues + self.lax_residues < len(border) or (self.lax_residues > self.confirmed_residues and self.confirmed_residues < (len(border) / 2)):  # If we have un-interacting atoms or the fitting was too bad. /how to go back and redo?/
                        cmake.superpose(chain, self.homo_chains, self.interactions_dict, self.chain_id_dict, self.new_pdb, border, self.tmp_count)
                        self.master.update_proteins([(self.new_pdb, 'tmp/part%s.pdb' % self.tmp_count)])
                    self.completed_borders += 1
                    self.interactions_checked = False

                else:
                    self.confirmed_residues, self.lax_residues = cmake.check_interactions(chain, border)
                    self.interactions_checked = True

                if self.completed_borders == len(self.interactions_dict[self.chain_id_dict[chain.get_id()]]):
                    self.completed_chains += 1
                    self.completed_borders = 0
                if self.completed_chains == len(self.new_pdb):
                    self.master.running = False
                    self.master.restart = True
                    self.master.run_butt['text'] = "Run "
                    self.master.log_label['text'] = "Job Finished Successfully"
                    self.master.output = self.new_pdb
                    sys.stdout.write("Job Finished Successfully\n")
                else:
                    self.run()
            else:
                while not self.master.running:
                    # to allow the function to be stopped but not killed we have this loop
                    # The sleep is there to reduce computational cost
                    # when the run button is pressed again this loop ends and the function calls itself again
                    time.sleep(1)
                self.run()


class PymolWindow(tk.Frame):
    def __init__(self, parent, master):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.master = master
        self.img = ImageTk.PhotoImage(Image.open("virus.jpg").resize((350, 350), Image.ANTIALIAS))
        self.panel = tk.Label(self, image=self.img)
        self.panel.pack()
        self.running = False
        self.queue = None
        self.job = None

    def update_image(self):
        self.img = ImageTk.PhotoImage(Image.open("tmp/tmp_img.png").resize((350, 350), Image.ANTIALIAS))
        self.panel['image'] = self.img

    def start_pymol_image_thread(self, pdb):
        if not self.running:
            self.running = True
            self.queue = queue.Queue()
            self.job = PymolImageThread(self.queue, pdb)
            self.job.start()
            self.master.current_module.log_label["text"] = "Making Picture"
            self.master.after(100, self.process_queue)

    def process_queue(self):
        try:
            msg = self.queue.get(0)
            if msg == "Pymol task finished\n":
                sys.stdout.write(msg)
                self.master.current_module.log_label["text"] = "Image processed Successfully"
                self.update_image()
            else:
                sys.stderr.write("Pymol not installed properly\n")
                self.master.current_module.log_label["text"] = "Error while processing Image"
            self.running = False
        except queue.Empty:
            self.master.after(100, self.process_queue)


class PymolImageThread(threading.Thread):
    def __init__(self, queue, pdb):
        threading.Thread.__init__(self)
        self.queue = queue
        self.pdb = pdb

    def run(self):
        """
        The reason why we launch a python script as a subprocess is because the cmd of pymol as a really strange behaviour
        for some reason it kills the main thread when it finish and also it does't wait until the process is finished before
        going to the next line.
        """
        result = subprocess.run(["python", "pymolmanager.py", "%s" % self.pdb])
        if result.returncode == 0:
            self.queue.put("Pymol task finished\n")
        else:
            self.queue.put("Pymol not installed properly\n")


class OutputOptions(tk.Frame):
    def __init__(self, parent, master):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.master = master
        self.info_butt = tk.Button(self, text='Info', command=lambda: self.master.console.change_channel('info'))
        self.stdout_butt = tk.Button(self, text='Stdout', command=lambda: self.master.console.change_channel('stdout'))
        self.stderr_butt = tk.Button(self, text='Stderr', command=lambda: self.master.console.change_channel('stderr'))
        self.clear_butt = tk.Button(self, text='Clear', command=self.master.console.clear_text)
        self.save_butt = tk.Button(self, text='Save', command=self.master.console.save)

        # Widgets placement after declaration
        self.info_butt.grid(row=0, column=0)
        self.stdout_butt.grid(row=0, column=1)
        self.stderr_butt.grid(row=0, column=2)
        self.clear_butt.grid(row=0, column=3)
        self.save_butt.grid(row=0, column=4)


class TextConsole(tk.Text):
    def __init__(self, parent, chanel, **kargs):
        tk.Text.__init__(self, parent, **kargs)
        self.parent = parent
        self.chanel = None
        self.text = {'stderr': ["Standard error:\n"],
                     'stdout': ["Standard Output:\n"],
                     'info': ["Protein Info:\n"]}
        self.change_channel(chanel)

    def clear_text(self):
        self.config(state='normal')
        self.delete(2.0, tk.END)
        self.config(state='disabled')
        self.text[self.chanel] = [self.text[self.chanel][0]]

    def update_text(self, chanel, data):
        self.text[chanel].append(data)
        if self.chanel == chanel:
            self.config(state='normal')
            self.insert(tk.END, data)
            self.config(state='disabled')

    def change_channel(self, chanel):
        if self.chanel != chanel:
            self.chanel = chanel
            self.config(state='normal')
            self.delete(1.0, tk.END)
            for line in self.text[self.chanel]:
                self.insert(tk.END, line)
            self.config(state='disabled')

    def save(self):
        out = filedialog.asksaveasfile(title="Save Stdout as:", mode='w', defaultextension=".txt")
        if out is not None:
            text2save = "".join(self.text['stdout'])
            out.write(text2save)
            out.close()
        error = filedialog.asksaveasfile(title="Save Stderr as:", mode='w', defaultextension=".txt")
        if error is not None:
            text2save = "".join(self.text['stderr'])
            error.write(text2save)
            error.close()


class PDBSplit(tk.Frame):
    def __init__(self, parent, master):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.master = master

class MainW(tk.Tk):
    def __init__(self, parent):
        # Basic Configuration
        tk.Tk.__init__(self, parent)
        self.parent = parent
        self.title("InteractMaker")
        self.current_module_name = "Interact Maker"
        # Menu
        self.menu = MenuBar(self)
        self.config(menu=self.menu)
        # Container Frames for the Widgets
        self.output_options_frame = tk.Frame(self, self.parent)
        self.output_frame = tk.Frame(self, self.parent)
        self.interact_frame = tk.Frame(self, self.parent, padx=5, width=807, height=355)
        self.pdbsplit_frame = tk.Frame(self, self.parent, width=807, height=355)
        self.pdbvalidate_frame = tk.Frame(self, self.parent, bg="blue", width=807, height=355)
        # Container Frames Placement
        self.interact_frame.pack(fill=tk.X)
        self.output_options_frame.pack(fill=tk.X)
        self.output_frame.pack(fill='both', expand=True)
        # Widgets Declaration
        self.current_module = MatchModule(self.interact_frame, self)
        self.img_display = PymolWindow(self.interact_frame, self)
        self.pdbsplitter = PDBSplit(self.pdbsplit_frame, self)
        self.console = TextConsole(self.output_frame, 'stderr', state='disabled')
        self.output_options = OutputOptions(self.output_options_frame, self)
        # Widgets Placement
        self.current_module.pack(side=tk.LEFT, anchor=tk.N)
        self.img_display.pack(side=tk.RIGHT)
        self.output_options.pack(fill='x')
        self.console.pack(fill='both', expand=True)


if __name__ == "__main__":
    if not os.path.exists('tmp'):
        os.mkdir('tmp')
    app = MainW(None)
    sys.stderr = Writer(app, 'stderr')
    sys.stdout = Writer(app, 'stdout')
    info_console = Writer(app, "info")
    from PDB import ProteinStructure as PS  # Imported here because we have overrode the sys.stderr and stdout channel
    import PDBaligner_up as cmake
    import temporal_cleaner as clean
    app.mainloop()
    clean.clean_temporal()
