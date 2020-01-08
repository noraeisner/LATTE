import tkinter as tk
#import tkSimpleDialog
from tkinter import simpledialog

ROOT = tk.Tk()
ROOT.withdraw()


# load first prompt window which asks for TIC ID    
#tic = simpledialog.askstring(title="TIC",
#                                      prompt="Enter TIC ID:")

class MyDialog(simpledialog.Dialog):

    def body(self, master):

        tk.Label(master, text="TIC ID:").grid(row=0)
        self.e1 = tk.Entry(master)
        self.e1.grid(row=0, column=1)
        
        #instructions = tk.Label(master, text="Check for FFI mode").grid(row=0)
        self.FFIbox = tk.IntVar()
        self.answer = tk.Checkbutton(master,text="Check for FFI mode", variable=self.FFIbox)
        self.answer.grid(row=1, column=1,  columnspan=2)

    def apply(self):
        global TIC
        global FFI

        ROOT.form=(self.FFIbox.get())
        TIC = (self.e1.get())
        FFI =  (self.FFIbox.get())
        

k = MyDialog(ROOT)





#class MyDialog(simpledialog.Dialog):
#
#    def body(self, master):
#
#        tk.Label(master, text="First:").grid(row=0)
#
#        self.e1 = tk.Entry(master)
#
#        self.e1.grid(row=0, column=1)
#
#        self.cb = tk.Checkbutton(master, text="Hardcopy")
#        self.cb.grid(row=2, columnspan=2)
#
#        return self.e1, self.cb # initial focus
#
#    def apply(self):
#        first = int(self.e1.get())
#
#k = MyDialog(ROOT)

