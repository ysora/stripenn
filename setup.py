import pathlib
import setuptools

HERE = pathlib.Path(__file__).parent

README = (HERE / "README.md").read_text()

setuptools.setup(
    name = 'stripenn',
    version = '1.0.7',
    author = 'Sora Yoon',
    author_email = 'sora.yoon@pennmedicine.upenn.edu',
    description = "Image-processing based detection of architectural stripes from chromatic conformation data",
    long_description=README,
    long_description_content_type="text/markdown",
    url = 'https://github.com/ysora/stripenn',
    packages_dir={'':'src'}
    packages=setuptools.find_packages('src'),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["cooler","pandas","numpy","matplotlib","opencv-python","scikit-image","scipy","joblib","tqdm","Bottleneck"],
    entry_points={
    	"console_script": ["stripenn=stripenn.stripenn:main"]
    },
)
