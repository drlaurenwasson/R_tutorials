<b> Intro to R Studio </b>
When you open R studio, you should see four windows:
- Top left window is source code. This is where your scripts will be shown to you.
- Top right is your local R environment. 
- Bottom left is the console where the actual code is run
- Bottom right is the folder where your files will be saved.

Start a new project using File -> New Project. Here you can specify the location where you want to run your project, and that will become your 
working directory. That directory will be shown in the bottom right (the Source file location). The working directory can be changed by going to "Session -> Change working directory

On the top left you should see a blank script titled Untitled1. Here is where you write your script

Say we wanted to print "TBX5" just on the screen. That would be

`print("TBX5")`

If you type that in the top left, and then hit Ctrl+Enter (Cmd+Enter on Mac), you will see "TBX5" show up on the Console (bottom left)

If you want to save "TBX5" as a variable to be used later, you would type:

`favgene<- "TBX5"`

This assigns "TBX5" to a variable called "favgene". If you then type 

`favgene` 

and run that, that will return "TBX5"

If you want to save more than one gene (say, a list of favorite genes):

`favgenes<- c("TBX5", "CHD4", "SMYD1")`

This will save your list of favorite genes as a list in your environment (top right) called "favgenes"

If you wanted to save this list as a table:

`write.table(favgenes, file = "favegenes.txt", quote = F)`

This will save your list as a file called favegenes.txt and will be stored in the directory shown on the bottom right.

