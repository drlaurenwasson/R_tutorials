# Introduction to R Studio

When you open R studio, you should see four windows:
- The window on the top left is your source code. This is where your scripts are located. 
- The window on the top right is your local R environment.
- The window on the bottom left is the console. This is where your actual code is run.
- The window on the bottom right should represent your working directory. This is where your files are saved.

Start a new project using File -> New Project. Here you can specify the location where you want to run your project, and that will become your working directory. That directory will be shown in the bottom right (the Source file location). The working directory can be changed by going to "Session -> Change working directory

On the top left you should see a blank script titled Untitled1. Here is where you write your script

Say we wanted to print "TBX5" just on the screen. That would be:


```R
print("TBX5")
```

    [1] "TBX5"
    

If you type that in the top left (source code), and then hit Ctrl+Enter (Cmd+Enter on Mac), you will see "TBX5" show up on the Console (bottom left).

If you want to save "TBX5" as a variable to be used later, you would type:


```R
favegene <- "TBX5"
```

This assigns "TBX5" to a variable called "favgene". If you then type:


```R
favegene
```
you will output "TBX5".
```
'TBX5'
```

If you want to save more than one gene (say, a list of favorite genes), you use the concatenate function: c()


```R
favegenes <- c("TBX5", "CHD4", "SMYD1")
```

This will save your list of favorite genes as a list in your environment (top right) called "favgenes"

If you wanted to save this list as a table:


```R
write.table(favegenes, file = "favegenes.txt", quote = F)
```

This will save your list as a file called favegenes.txt and will be stored in the directory shown on the bottom right.
