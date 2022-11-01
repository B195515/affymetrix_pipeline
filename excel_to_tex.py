import pandas as pd
filename = input("Enter filename: ")
df = pd.read_excel(f'{filename}.xls')
df.to_latex(f'{filename}.tex')
print(f'{filename}.tex is created')
