# Handy R tricks and tips

These are some R tips that are not specific to any package, but more for R in general.

## Manipulating data frames

The data frame is likely the most useful thing about R, it can handle large data sets very easily. 


### Rows and columns

Say we have a data frame of an annotated bed file called ```peaks.db.PGP1enriched```. 
First, to view the data frame in a separate R window: 

```
View(peaks.db.PGP1enriched)
```

To see row and column names
```
row.names(peaks.db.PGP1enriched)
colnames(peaks.db.PGP1enriched)
```

#### Isolating rows and columns 

To take all the rows, but the first three columns of a data frame and assign it to a new variable:
```
newdf<- peaks.db.PGP1enriched[,1:3]
```
To take all the columns, but the first three rows of a data frame and assign it to a new variable:

```
newdf<- peaks.db.PGP1enriched[1:3,]
```

The rows are entered before the comma, and the columns are entered after. If you wanted just the first three rows and columns 2-5, it would look like this

```
newdf<- peaks.db.PGP1enriched[1:3,2:5]
```

Lastly, if you want all rows, but columns 1,3,5,and 8, it would look like this:
```
newdf<- peaks.db.PGP1enriched[,c(1,3,5,8)]
```

### Subsetting by certain values in a dataframe:
Coming soon!
