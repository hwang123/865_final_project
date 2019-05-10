import os
import sys
import random
import math
import numpy as np
import skimage.io
import matplotlib
import matplotlib.pyplot as plt

# Root directory of the project
ROOT_DIR = os.path.abspath("../")

# Import Mask RCNN
sys.path.append(ROOT_DIR)  # To find local version of the library
from mrcnn import utils
import mrcnn.model as modellib
from mrcnn import visualize
# Import COCO config
sys.path.append(os.path.join(ROOT_DIR, "samples/coco/"))  # To find local version
import coco

# %matplotlib inline 

def write_results_to_file(r): 
	# Write out the results of MaskRCNN 
	rois = r['rois']
	class_ids = r['class_ids']
	scores = r['scores']
	masks = r['masks']

	# Image dimensions
	height = len(masks)
	width = len(masks[0])
	# Number of ROIs
	num_regions = len(rois)

	# Open file
	f = open("../../res/coco_res.txt","w+")

	# Write out the width and height of the image on one line
	f.write("%d " % width)
	f.write("%d\n" % height)

	# Write out the number of ROI on one line
	f.write("%d\n" % (num_regions))

	for i in range(num_regions):
		class_id = class_ids[i];
		roi = rois[i]

		score = scores[i]

		# Write out the class ID and prediction score on one line
		f.write("%d " % class_id)
		f.write("%f\n" % score)

		# Write out the ROI bounding box on one line
		for j in range(3):
			coord = roi[j]
			f.write("%d " % coord)
		f.write("%d\n" % roi[3])

		# Write out the pixel mask of the ROI
		x1 = roi[1]
		x2 = roi[3]

		y1 = roi[0]
		y2 = roi[2]

		for y in range(y1, y2 + 1):
			for x in range(x1, x2 + 1):
				mask_val = masks[y][x][i]
				pixel = 1 if mask_val else 0;
				f.write("%d " % pixel)
			f.write("\n")

	f.close()


# Directory to save logs and trained model
MODEL_DIR = os.path.join(ROOT_DIR, "logs")

# Local path to trained weights file
COCO_MODEL_PATH = os.path.join(ROOT_DIR, "mask_rcnn_coco.h5")
# Download COCO trained weights from Releases if needed
if not os.path.exists(COCO_MODEL_PATH):
    utils.download_trained_weights(COCO_MODEL_PATH)

# Directory of images to run detection on
IMAGE_DIR = os.path.join(ROOT_DIR, "images")

# Configurations
class InferenceConfig(coco.CocoConfig):
    # Set batch size to 1 since we'll be running inference on
    # one image at a time. Batch size = GPU_COUNT * IMAGES_PER_GPU
    GPU_COUNT = 1
    IMAGES_PER_GPU = 1

config = InferenceConfig()
# config.display()

# Create model object in inference mode.
model = modellib.MaskRCNN(mode="inference", model_dir=MODEL_DIR, config=config)

# Load weights trained on MS-COCO
model.load_weights(COCO_MODEL_PATH, by_name=True)

# COCO Class names
# Index of the class in the list is its ID. For example, to get ID of
# the teddy bear class, use: class_names.index('teddy bear')
class_names = ['BG', 'person', 'bicycle', 'car', 'motorcycle', 'airplane',
               'bus', 'train', 'truck', 'boat', 'traffic light',
               'fire hydrant', 'stop sign', 'parking meter', 'bench', 'bird',
               'cat', 'dog', 'horse', 'sheep', 'cow', 'elephant', 'bear',
               'zebra', 'giraffe', 'backpack', 'umbrella', 'handbag', 'tie',
               'suitcase', 'frisbee', 'skis', 'snowboard', 'sports ball',
               'kite', 'baseball bat', 'baseball glove', 'skateboard',
               'surfboard', 'tennis racket', 'bottle', 'wine glass', 'cup',
               'fork', 'knife', 'spoon', 'bowl', 'banana', 'apple',
               'sandwich', 'orange', 'broccoli', 'carrot', 'hot dog', 'pizza',
               'donut', 'cake', 'chair', 'couch', 'potted plant', 'bed',
               'dining table', 'toilet', 'tv', 'laptop', 'mouse', 'remote',
               'keyboard', 'cell phone', 'microwave', 'oven', 'toaster',
               'sink', 'refrigerator', 'book', 'clock', 'vase', 'scissors',
               'teddy bear', 'hair drier', 'toothbrush']

# Load a random image from the images folder
file_names = next(os.walk(IMAGE_DIR))[2]

image = skimage.io.imread(os.path.join(IMAGE_DIR, '2516944023_d00345997d_z.jpg'))

# Run detection
results = model.detect([image], verbose=1)
r = results[0]

# Write results to file
write_results_to_file(r)

# Visualize results
if False:
	visualize.display_instances(image, r['rois'], r['masks'], r['class_ids'], 
                            class_names, r['scores'])
