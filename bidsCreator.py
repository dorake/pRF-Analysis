"""
-----------------------------------------------------------------------------------------
bidsMBsequence.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create bids formated data out of seq_test experiment collected on TK
-----------------------------------------------------------------------------------------
Input(s):
none
-----------------------------------------------------------------------------------------
Output(s):
=> formated bids files
-----------------------------------------------------------------------------------------
To run:
python /home/shared/2018/visual/seq_test/codes/post_cart.py
-----------------------------------------------------------------------------------------
"""
import numpy as np
import json 
import os, re

def create_jsonFiles(b0IntendedForNames,subName):
	#sliceTiming_MB1 = np.linspace(0, 1.5, 16)
	sliceTiming_MB1 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
	#sliceTiming_MB1 = sliceTiming_MB_1[:-1].tolist()
	sliceTiming_MB2 = np.concatenate((sliceTiming_MB1, sliceTiming_MB1)).tolist()
	sliceTiming_MB3 = np.concatenate((sliceTiming_MB2, sliceTiming_MB1)).tolist()
	sliceTiming_MB4 = np.concatenate((sliceTiming_MB3, sliceTiming_MB1)).tolist()

	## Fill with experiment information, e.g.
	jsonTask_MB1Int128 = {
							"EchoTime": ,
							"EpiFactor": ,
							"SenseFactor": ,
							"PhaseEncodingDirection": "j",
							"SliceDirection": "z",
							"RepetitionTime": ,
							"MultiBandFactor": ,
							"NumberDummyScan": ,
							"WaterFatShift": ,
							"EchoTimeDifference": , 
							"EffectiveEchoSpacing": , 
							"SliceTiming": sliceTiming_MB1,
							"SliceOrder": "ascending", 
							"TaskName": "MB1Int128"
							}

	fmap_json = {
				"EchoTime1": ,
				"EchoTime2": ,
				"IntendedFor":b0IntendedForNames
				}

	jsonFileNames = [jsonTask_MB1Int128, jsonTask_MB1Int192, jsonTask_MB1Int256,\
					 jsonTask_MB2Int128, jsonTask_MB2Int192, jsonTask_MB2Int256,\
					 jsonTask_MB3Int128, jsonTask_MB3Int192, jsonTask_MB3Int256,\
					 jsonTask_MB4Int128, jsonTask_MB4Int192, jsonTask_MB4Int256]

	jsonPath = '/home/shared/2018/visual/seq_test/bids/'
	fmapPath = '/home/shared/2018/visual/seq_test/bids/'+ subName +'/fmap/'

	for file in jsonFileNames:
		taskName = file['TaskName']

		completeName = jsonPath + 'task-' + taskName + '_bold.json'
		
		with open(completeName, 'w') as outfile:
			json.dump(file, outfile)


	b0JsonName = fmapPath + subName + '_phasediff.json'

	with open(b0JsonName, 'w') as outfile:
		json.dump(fmap_json, outfile)

	return None

def copy_file(srcDir, destDir):

	if os.path.isfile(destDir) == False:
		
		try:
			working_string = 'rsync -az --no-g --no-p --progress ' + srcDir + ' ' + destDir
			#print(working_string) # use this line to debug
			os.system(working_string)
		except OSError:
			print("File already exists: " + destDir)
			print("")
	else:
		print("File already exists: " + destDir)

	return None


def create_dictionary(destDir):

	if os.path.isdir(destDir) == False:

		try:
			os.makedirs(destDir)
		except OSError:
			Print('Directory already exists:' + destDir)
			print("")
	else:
		print("Directory already exists: " + destDir)
	return None	

def getFilesFromDict(subjDir, regexArgument):

	## Get all files in the nifti directory saved in subjNiftis
	filesInDir = []
	for (dirpath, dirnames, filenames) in os.walk(subjDir):
		for file in filenames:
			filesInDir.append(file)
		break

	## Extract all relevant files from directory by using regex
	fileList = []
	for file in filesInDir:
		if re.search(regexArgument, file):
			fileList.append(file)

	#print(fileList)
	return fileList


original_data = ''
expFolder = ''
subName = 'sub-001'
subjFolder = expFolder + '/' + subName
fmapFolder = subjFolder + '/fmap'
funcFolder = subjFolder + '/func'
anatFolder = subjFolder + '/anat'

nrMBCorrespondence = {
						5: 'MB1',
						6: 'MB2',
						8: 'MB3',# note: run 7 was a mistake, bold signal data were deleted but not behavioral or eyeData
						10: 'MB4'
						}

edfMBCorrespondence = {
						'1': 5,
						'2': 6,
						'4': 8,	# note: run 7 was a mistake, bold signal data were deleted but not behavioral or eyeData
						'5': 10
						}

create_dictionary(expFolder)
create_dictionary(subjFolder)
create_dictionary(fmapFolder)
create_dictionary(funcFolder)
create_dictionary(anatFolder)


oldFuncFolder = ''
native_funcRegex = 'SENSE_([0-9]|[0-9][0-9])_1.nii.gz'
int192_Regex = 'Interp192_([0-9]|[0-9][0-9])_[0-9].nii.gz'
int256_Regex = 'Interp256_([0-9]|[0-9][0-9])_[0-9].nii.gz'
edf_Regex = 'tk_[0-9]_.*.edf'
pickle_Regex = 'tk_[0-9]_.*.pickle'

nativeFileNames = getFilesFromDict(oldFuncFolder,native_funcRegex)
int192FileNames = getFilesFromDict(oldFuncFolder, int192_Regex)
int256FileNames = getFilesFromDict(oldFuncFolder, int256_Regex)
edfFileNames = getFilesFromDict(oldFuncFolder, edf_Regex)
pickleFileNames = getFilesFromDict(oldFuncFolder, pickle_Regex)

fileList = []
b0IntendedForNames = []

for file in int192FileNames:
	oldFileName = oldFuncFolder + '/' + file

	splittedFileName = file.split('_')
	fileNr = int(splittedFileName[-2])

	fileName = funcFolder + '/' + subName + '_task-' + nrMBCorrespondence[fileNr] + 'Int192_bold.nii.gz'
	fileTuple = (oldFileName, fileName)
	b0FileName = 'func/'+ subName + '_task-' + nrMBCorrespondence[fileNr] + 'Int192_bold.nii.gz'
	b0IntendedForNames.append(b0FileName)

	fileList.append(fileTuple)


for file in int256FileNames:
	oldFileName = oldFuncFolder + '/' + file

	splittedFileName = file.split('_')
	fileNr = int(splittedFileName[-2])

	fileName = funcFolder + '/' + subName + '_task-' + nrMBCorrespondence[fileNr] + 'Int256_bold.nii.gz'
	fileTuple = (oldFileName, fileName)
	b0FileName = 'func/' + subName + '_task-' + nrMBCorrespondence[fileNr] + 'Int256_bold.nii.gz'
	b0IntendedForNames.append(b0FileName)
	
	fileList.append(fileTuple)


for file in nativeFileNames:
	oldFileName = oldFuncFolder + '/' + file

	splittedFileName = file.split('_')
	fileNr = int(splittedFileName[-2])

	fileName = funcFolder + '/' + subName + '_task-' + nrMBCorrespondence[fileNr] + 'Int128_bold.nii.gz'
	fileTuple = (oldFileName, fileName)
	b0FileName = 'func/' + subName + '_task-' + nrMBCorrespondence[fileNr] + 'Int128_bold.nii.gz'
	b0IntendedForNames.append(b0FileName)

	fileList.append(fileTuple)

for file in edfFileNames:
	oldFileName = oldFuncFolder + '/' + file

	splittedFileName = file.split('_')
	fileNr = str(splittedFileName[1])

	if fileNr != '3':
		edfMB = edfMBCorrespondence[fileNr]
		
		fileName1 = funcFolder + '/' + subName + '_task-' + nrMBCorrespondence[edfMB] + 'Int128_eyeData.edf'
		fileName2 = funcFolder + '/' + subName + '_task-' + nrMBCorrespondence[edfMB] + 'Int192_eyeData.edf'
		fileName3 = funcFolder + '/' + subName + '_task-' + nrMBCorrespondence[edfMB] + 'Int256_eyeData.edf'
		
		fileTuple1 = (oldFileName, fileName1)
		fileTuple2 = (oldFileName, fileName2)
		fileTuple3 = (oldFileName, fileName3)

		fileList.append(fileTuple1)
		fileList.append(fileTuple2)
		fileList.append(fileTuple3)
	else:
		pass

for file in pickleFileNames:
	oldFileName = oldFuncFolder + '/' + file

	splittedFileName = file.split('_')
	fileNr = str(splittedFileName[1])

	if fileNr != '3':
		edfMB = edfMBCorrespondence[fileNr]
		
		fileName1 = funcFolder + '/' + subName + '_task-' + nrMBCorrespondence[edfMB] + 'Int128_behav.pickle'
		fileName2 = funcFolder + '/' + subName + '_task-' + nrMBCorrespondence[edfMB] + 'Int192_behav.pickle'
		fileName3 = funcFolder + '/' + subName + '_task-' + nrMBCorrespondence[edfMB] + 'Int256_behav.pickle'
		
		fileTuple1 = (oldFileName, fileName1)
		fileTuple2 = (oldFileName, fileName2)
		fileTuple3 = (oldFileName, fileName3)

		fileList.append(fileTuple1)
		fileList.append(fileTuple2)
		fileList.append(fileTuple3)

## Copying anatomical file
oldAnatFileName = '/home/fs_subjects/TK_290615/mri/T1.nii.gz'
newAnatFileName = anatFolder + '/' + subName + '_T1w.nii.gz'
anatTuple = (oldAnatFileName, newAnatFileName)
fileList.append(anatTuple)

## Copying field map files
tk_b0_magnitude = original_data + 'TK_B0/MB_MS_R6_WIP_B0_magnitude1.nii.gz'
tk_b0_phasediff = original_data + 'TK_B0/MB_MS_R6_WIP_B0_phasediff.nii.gz'

newB0MagnitudeName = fmapFolder + '/' + subName + '_magnitude1.nii.gz'
newB0phaseDiffName = fmapFolder + '/' + subName + '_phasediff.nii.gz'

magnitudeTuple = (tk_b0_magnitude, newB0MagnitudeName)
phaseDiffTuple = (tk_b0_phasediff, newB0phaseDiffName)
fileList.append(magnitudeTuple)
fileList.append(phaseDiffTuple)

## Create JSON file with dataset information
dataSetDescription = {
						'Name' : 'sequence_test',
						'BIDSVersion' : '1.0.2'
						}

datasetJson = expFolder + '/' + 'dataset_description.json'
with open(datasetJson, 'w') as outfile:
	json.dump(dataSetDescription, outfile)

create_jsonFiles(b0IntendedForNames,subName)

## Doing the actual file copying of all created tuple pairs (oldFile, newFile)
for tpl in fileList:
	copy_file(tpl[0], tpl[1])