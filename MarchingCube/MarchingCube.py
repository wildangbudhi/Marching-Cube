from MarchingCubeCpp import MarchingCube, ConvertHu
import numpy as np
from pydicom import read_file
from os import listdir
from time import time

def Make3D(DicomSeriesPath: str, DicomType: str, Threshold: int = 400):

	print("Load")
	slices = [read_file(DicomSeriesPath + file) for file in listdir(DicomSeriesPath)]
	slices.sort(key = lambda x: float(x.ImagePositionPatient[2]))

	pixel_Array = np.array([s.pixel_array for s in slices], dtype=np.float)
	ZDist = slices[1].SliceLocation - slices[0].SliceLocation
	YDist = slices[0].PixelSpacing[0]
	XDist = slices[0].PixelSpacing[1]

	if(DicomType == "CT Scan"):
		print("ConvertHu Start")
		intercept = pixel_Array[0].RescaleIntercept if hasattr(pixel_Array[0], 'RescaleIntercept') else -1024
		slope = pixel_Array[0].RescaleSlope if hasattr(pixel_Array[0], 'RescaleSlope') else 1
		ConvertHu(pixel_Array, slope, intercept)

	start = time()
	MarchingCube(pixel_Array, pixel_Array.shape[0], pixel_Array.shape[1], pixel_Array.shape[2], ZDist, YDist, XDist, Threshold, 2, DicomType + ".obj", True)
	print("Sukses:", time() - start)

def main():
	Make3D("D:/Python/3D Brain/Editing Code/Azis/Azis CT Scan/A/", "CT Scan", Threshold=400)
	#Make3D("D:/Python/3D Brain/Editing Code/Azis/Azis MRI/A/", "MRI", Threshold=400)
	Make3D("D:/Python/3D Brain/Editing Code/Azis/Azis MRA & DTI/501/", "MRA", Threshold=1000)

if __name__ == "__main__":
	main()