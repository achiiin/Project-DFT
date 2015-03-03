"""
Created on Tue Mar 03 15:57:41 2015
@author: Ching-Yen
"""
"""
This script is to get the Z-Score for each building block of SMILES of outliers.

input:
argument: percentage of outliers

 
""" 
### the HelpFormatter and provide a special intro for the options that should be handled "raw"
### Any other calls to .add_argument() where the help does not start with R| will be wrapped as normal.
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        # this is the RawTextHelpFormatter._split_lines
        if text.startswith('R|'):
            return text[2:].splitlines()  
        return argparse.HelpFormatter._split_lines(self, text, width)
        

       
parser = argparse.ArgumentParser(description="Do Linear Regression to the data", 
                                 formatter_class=SmartFormatter)
                           
                                 
                              
parser.add_argument("property", type=int, choices=[1,2,3,4],
                    help="R|choose the options of property:\n"
                    "1 = HOMO  2 = LUMO  3 = Dipole  4 = Gap")
parser.add_argument("percent_outlier", type=float, help="R|define the outlier\n"
                    "enter the percentage of data taken out to be outliers\n"
                    "input only number part")
parser.add_argument("-x", "--x_flavor", type=int, 
                    choices=[1,2], help="R|choose the flavor for x-axis:\n"
                    "1 = BP86/SVP    2 = B3LYP/SVP")
parser.add_argument("-y", "--y_flavor", type=int, 
                    choices=[1,2,3,4,5,6,7,8,9], help="R|choose the flavor for y-axis:\n"
                    "1 = BP86/SVP    2 = B3LYP/SVP  3 = PBE0/SVP \n4 = BH&HLYP/SVP 5 = M06-2X/SVP"
                    " 6 = HF/SVP \n7 = BP86/TZVP   8 = B3LYP/TZVP 9 = PBE0/TZVP")
                    
parser.add_argument("-rb","--result_box",action='store_false',default=True,
                    help="Hide the result display box in all plots")
parser.add_argument("-rb1","--result_box1",action='store_false',default=True,
                    help="Hide the result display box only in original data pool")
parser.add_argument("-rb2","--result_box2",action='store_false',default=True,
                    help="Hide the result display box only in outlier")
parser.add_argument("-rb3","--result_box3",action='store_false',default=True,
                    help="Hide the result display box only in remaining pool")
                    
parser.add_argument("-lr","--lr_line",action='store_false',default=True,
                    help="Hide all the LinearReg lines.")
parser.add_argument("-lr1","--lr_line1",action='store_false',default=True,
                    help="Hide the LinearReg line for original data pool.")
parser.add_argument("-lr2","--lr_line2",action='store_false',default=True,
                    help="Hide the LinearReg line for outliers.")
parser.add_argument("-lr3","--lr_line3",action='store_false',default=True,
                    help="Hide the LinearReg line for remaining pool.")
                   
args = parser.parse_args('1 1 -x 1 -y 6 -lr2 -lr3'.split())
#args = parser.parse_args()