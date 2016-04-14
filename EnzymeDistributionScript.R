##############################################################################
###
### Graphic visualisation of fluorescence absorbance
### Microplate reader: Infinite M200 (TECAN, Austria) 
###
############################################################################## 

################################################
## global variables
ReadingError=FALSE
TimeForSet1=0
TimeForSet2=0


################################################
## conversion of negative values to positive

NonNegative = function(Matrix)
  {
  Matrix[Matrix<0] = 0
  return(Matrix)
  }


################################################
## test the file suffix *.csv

GetSuffix = function(fileName) 
	{
	correct=TRUE 
	ChNum = nchar(fileName)
	s=strsplit(fileName, NULL)[[1]]
	a=ChNum-2
	suffix=s[a:ChNum]
	suffix=paste(c(suffix[1],suffix[2],suffix[3]), collapse="") 
	csv=strsplit('csv', NULL)[[1]]	
	csv=paste(c(csv[1],csv[2],csv[3]), collapse="") 
	if (csv!=suffix) correct=FALSE
	
	if (correct==FALSE) 
		{
		info = paste('Wrong file suffix - file ignored > ',fileName, collapse="") 
		print(info)
		}
	return(correct)
	}

################################################
## load file

FileLoad = function(fileName) 
	{
	result=try(read.delim(fileName,sep=',')) 
	if (inherits(result, 'try-error')) 
		{
		print('Chyba ve cteni souboru!!!')
		ReadingError=TRUE
		}
		else
		{
		mat=read.delim(fileName,sep=';')	
		info=c('nacteno: ',fileName)
		print(info)
		}
	return(mat)
	}


################################################
## interpolation
# value of absorbance for specific time 
# based on two consecutive measurement

TimeInterpolation = function(StartTimeM1, M1, StartTimeM2, M2, Interval)
  {
  maxCol=ncol(M1)
  maxRow=nrow(M1)
  
  # tvorba nove matice	
  InterPolMat = matrix(ncol = maxCol, nrow = maxRow)
  
  # spocte hodnotu kdy byl zmeren posledni vzorek
  Tmax= StartTimeM1 + maxCol*maxRow*Interval

  x=1
  y=1

  poz=1
  back=FALSE	# direction of matrix walktrough 
                # data are alternately reading from left to right and contrariwise

  while (y <= maxRow)
      {

      if (back==FALSE)
        {
        x=1
        while (x <= maxCol)
          {
          T1=StartTimeM1+poz*Interval
          T2=StartTimeM2+poz*Interval

          Proportion=(Tmax-T1)/(T2-T1)
          InterPolMat[y,x]= M1[y,x]+ (Proportion*(M2[y,x]-M1[y,x]))

          x = x+1
          poz = poz+1
          back=TRUE
          }
        }
        else
        {
        x=maxCol
        while (x > 0)
          {
          T1=StartTimeM1+poz*Interval
          T2=StartTimeM2+poz*Interval

          Proportion=(Tmax-T1)/(T2-T1)
          InterPolMat[y,x]= M1[y,x]+ (Proportion*(M2[y,x]-M1[y,x]))

          x = x-1
          poz = poz+1
          back=FALSE
          }
        }

      y = y+1
      }
  return(InterPolMat)
  }


################################################
## draw

DrawDiference = function(Mat1,Mat2,FilesNames,OutPutPath,order,SaveResult,ExperimentName)
{
#print(Mat1)
#print(Mat2)
DifMat = NonNegative(Mat2-Mat1)
#print(DifMat)

savingCSVFile=paste(OutPutPath,'/output',order,'.csv',sep="") 
write.table(DifMat,file=savingCSVFile,sep=",",
            quote=FALSE,row.names=FALSE,col.names=FALSE)


x = 10*(1:nrow(DifMat))
y = 10*(1:ncol(DifMat))  


maxA = max(DifMat)
minA = min(DifMat)
progress = round((max(DifMat)-min(DifMat))/15)

savingFile=paste(OutPutPath,'/output',order,'.png',sep="") 


#image(x, y, DifMat, col = terrain.colors(100), axes = FALSE)


image(x, y, DifMat, col = terrain.colors(100), axes = FALSE)
box()



if (SaveResult) {
		png(filename=savingFile, width = 1000, height = 1000, bg = "white", pointsize = 22) 
		print(paste('ukladam...',savingFile))
		}



filled.contour(x, y, DifMat, 
	
    color = terrain.colors,
    plot.title = title(main = ExperimentName, sub = '',
      cex.main = 1.5,   font.main= 1, col.main= "black",
      cex.sub = 1.5, font.sub = 1, col.sub = "blue"),
    key.title = title(main = list("Absorbance\n(relative)", cex=0.90,
                  col="black", font=1)),
plot.axes = { box() },


    key.axes = axis(4, seq(minA, maxA, by = progress)),asp=1,axes = FALSE) 


# info
mtext(FilesNames,side = 1, line = 2,  cex = .99)
mtext(paste("T. Vetrovsky 2010 (c) - ", R.version.string),
      	side = 1, line = 4, adj = 1, cex = .66)




if (SaveResult) dev.off()
}



################################################
## Main Function
# 
# InPutPath 	      - input data path (default = D:/)
# MeasurementInterval - time difference between each measurements
# ScanInterval        - time of single value measurement
# OutPutPath          - output data path

Main = function(ExperimentName, InPutPath, MeasurementInterval, ScanInterval, OutPutPath)
{
if (InPutPath=='default')
	{
	InPutPath='D:/'
	setwd(InPutPath)
	}
	else
	{setwd(InPutPath)}

# save the resulted graph
SaveResult=TRUE
if (OutPutPath=='no')SaveResult=FALSE
if (OutPutPath=='default') OutPutPath=InPutPath

# zjisti pocet souboru v cilove slozce
files = list.files(path = ".", pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE)

MatNumber=length(files)

# ST - pocatecni cas, MI - rozdil mezi jednotlivymi merenimi
ST=0
MI=MeasurementInterval

if (MatNumber>3) 
	{ 
	## predprivaneni "listu" matic
	MatrixList=vector("list",MatNumber)
	FilesNamesList=vector("list",MatNumber)
	Count=0
	for (x in 1:MatNumber) 
		{	
		if (GetSuffix(files[x]))
			{
			FilesNamesList[[Count+1]] = files[x]
			MatrixList[[Count+1]]=as.matrix(NonNegative(FileLoad(files[x]))) 
			Count=Count+1
			print(paste('ok',files[x]))
			}
			
		}

	print(Count)

	if (ReadingError!=TRUE)
		{
		if (Count>3){
		for (x in 1:(Count-3)) 
			{	
			
			T1=MatrixList[[x]]
			T2=MatrixList[[x+1]]
			T3=MatrixList[[x+2]]
			T4=MatrixList[[x+3]]


			IPolMat1 = TimeInterpolation(ST+(x*MI), T1, ST+((x+1)*MI), T2,ScanInterval) 
			IPolMat2 = TimeInterpolation(ST+((x+2)*MI), T3, ST+((x+3)*MI), T4, ScanInterval)
			FilesNames = paste('uses files: ',FilesNamesList[[x]],',',FilesNamesList[[x+1]],',',
			FilesNamesList[[x+2]],',',FilesNamesList[[x+3]],collapse="") 
			DrawDiference(IPolMat1,IPolMat2,FilesNames,OutPutPath,x,SaveResult,ExperimentName)
			}

	
		print('finished')
		}
		else print('too few files for analysis...minimum is 4 *.csv files')
		}
		else print('some error occurred!!!')
	
	}
	else
	print('too few files for analysis...minimum is 4 *.csv files')
}


################################################
## RUN SCRIPT EXAMPLE

Main("CHITINASE",'D:/data/input',600,0.3645833,'D:/data/output')
