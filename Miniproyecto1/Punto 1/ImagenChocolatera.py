import cv2
import numpy as np
import matplotlib.pyplot as plt

# Load the image with alpha channel
image = cv2.imread('Chocolatera.png', cv2.IMREAD_UNCHANGED)  # This gives BGRA

# Split the channels
b, g, r, a = cv2.split(image)

# Convert BGRA to RGBA for matplotlib
rgba = cv2.merge([r, g, b, a])

# Create a white background
white_bg = np.ones_like(rgba[:, :, :3], dtype=np.uint8) * 255

# Normalize alpha channel to [0,1]
alpha = a.astype(float) / 255.0
alpha = np.repeat(alpha[:, :, np.newaxis], 3, axis=2)

# Alpha blending
blended = (rgba[:, :, :3] * alpha + white_bg * (1 - alpha)).astype(np.uint8)

x = [929, 920, 908, 893, 885, 881, 880, 879, 880, 887, 912, 917, 936, 953, 943, 929, 891,
     877, 863, 844, 823, 799, 781, 757, 727, 699, 671, 634, 948, 912, 951, 900, 952, 899]
y = [94, 110, 135, 178, 213, 240, 250, 269, 327, 370, 449, 462, 511, 645, 740, 786, 856,
     874, 890, 909, 926, 941, 950, 959, 967, 972, 975, 976, 566, 822, 688, 415, 607, 412]
x = np.array(x)
y = np.array(y)
x_corrected = x - 634
y_corrected = y - 94

y_sorted = np.sort(y_corrected)
x_sorted = np.array([0]*len(y_sorted))
order = [0]*len(y_sorted)

for i in range(len(y_sorted)):
    order[i] = int(np.where(y_corrected == y_sorted[i])[0][0])
    x_sorted[i] = x_corrected[order[i]]

y_scaled = np.array(x_sorted*0.2081218274)
x_scaled = np.array(y_sorted*0.2081218274)
print("x = ", x_scaled)
print("y = ", y_scaled)

# Display
plt.imshow(blended)
plt.scatter(x, y)
plt.axis('off')
plt.title('Blended PNG on White Background')
plt.show()

plt.scatter(x_corrected, y_corrected)
plt.scatter(x_scaled, y_scaled)
plt.axis()
plt.grid()
plt.show()
