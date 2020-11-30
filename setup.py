import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

README = (HERE / "README.md").read_text()

setup(
    name = 'stripenn',
    version = '1.1.5',
    author = 'Sora Yoon',
    author_email = 'sora.yoon@pennmedicine.upenn.edu',
    description = "Image-processing based detection of architectural stripes from chromatic conformation data",
    long_description=README,
    long_description_content_type="text/markdown",
    url = 'https://github.com/ysora/stripenn',
    #package_dir={'': 'src'},
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["cooler","pandas","numpy","matplotlib","opencv-python","scikit-image","scipy","joblib","tqdm",'typer','pathlib'],
    entry_points={
    	"console_scripts": ["stripenn=stripenn.cli:main"]
    },
)
