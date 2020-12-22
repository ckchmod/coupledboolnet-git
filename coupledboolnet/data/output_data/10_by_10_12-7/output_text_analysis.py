import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv("output.txt")
df = df.rename(columns={"output.txt":"output"})

df = df.output.str.split("_", expand=True)
df[6] = df[6].str.rstrip(".pkl")

newcols  = {0:"sys_arg", 1:"indexcount", 2:"simcount", 3:"k", 4:"p", 5:"h", 6:"T_c"}

df = df.rename(columns=newcols)

cols = ["k", "p", "h", "T_c"]
df[cols] = df[cols].apply(pd.to_numeric)
df["lyapunov"] = 2*df["k"]*df["p"]*(1-df["p"])
print(df)

df.plot.scatter(x="lyapunov", y="h")
plt.show()
