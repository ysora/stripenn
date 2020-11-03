import pathlib
import setuptools

HERE = pathlib.Path(__file__).parent

README = (HERE / "README.md").read_text()

setuptools.setup(
    name = 'stripenn',
    version = '1.0.0',
    description = "Image-processing based detection of architectural stripes from chromatic conformation data",
    long_description=README,
    long_description_content_type="text/markdown",
    scripts = ['stripenn'],
    author = 'Sora Yoon',
    author_email = 'sora.yoon@pennmedicine.upenn.edu',
    description = 'Architectural stripe detection program',
    long_description = "Architectural stripe detection tool using image processing methods",
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/VahediLab/stripenn',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=['stripenn']
    include_package_data=False,
    install_requires=["argparse","cooler","multiprocessing","pandas",
    	"numpy","math","matplotlib","opencv-python","statistics",
    	"scikit-image","scipy","joblib","tqdm","Bottleneck"],
    entry_points={
    	"console_script": [
    		"stripenn=src.stripenn:main"
    	]
    },
)
