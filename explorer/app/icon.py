import os
import subprocess
from PIL import Image

def create_iconset(png_path, iconset_dir):
    """
    Creates an iconset folder with various sizes from a single PNG image.
    """
    if not os.path.exists(iconset_dir):
        os.makedirs(iconset_dir)
    
    # Define the sizes you want to generate.
    # Common sizes for macOS icons include 16, 32, 128, 256, and 512.
    sizes = [
        (16, 16),
        (32, 32),
        (128, 128),
        (256, 256),
        (512, 512)
    ]
    
    # Open the original PNG image.
    base_image = Image.open(png_path)
    
    # For each size, create a normal image and a retina (@2x) version.
    for w, h in sizes:
        # Normal image
        normal_size = (w, h)
        normal_img = base_image.resize(normal_size, Image.LANCZOS)
        normal_filename = f"icon_{w}x{h}.png"
        normal_img.save(os.path.join(iconset_dir, normal_filename))
        
        # Retina image (double the size)
        retina_size = (w * 2, h * 2)
        retina_img = base_image.resize(retina_size, Image.LANCZOS)
        retina_filename = f"icon_{w}x{h}@2x.png"
        retina_img.save(os.path.join(iconset_dir, retina_filename))
    
    print(f"Created iconset folder '{iconset_dir}' with resized images.")

def convert_iconset_to_icns(iconset_dir, icns_path):
    """
    Uses the macOS 'iconutil' command to convert an iconset folder into an ICNS file.
    """
    try:
        subprocess.run(["iconutil", "-c", "icns", iconset_dir, "-o", icns_path],
                       check=True)
        print(f"ICNS file created at: {icns_path}")
    except subprocess.CalledProcessError as e:
        print("Error converting iconset to ICNS:", e)

if __name__ == '__main__':
    # Input PNG file (adjust the filename/path as needed)
    png_file = "icn.png"
    
    # Folder where the iconset will be created
    iconset_folder = "icon.iconset"
    
    # Output ICNS file
    icns_file = "icn.icns"
    
    # Create the iconset folder with all required sizes
    create_iconset(png_file, iconset_folder)
    
    # Convert the iconset folder to an ICNS file using iconutil
    convert_iconset_to_icns(iconset_folder, icns_file)