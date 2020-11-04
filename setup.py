import pathlib
import setuptools

HERE = pathlib.Path(__file__).parent

README = (HERE / "README.md").read_text()

setuptools.setup(
    name = 'stripenn',
    version = '1.0.4',
    description = "Image-processing based detection of architectural stripes from chromatic conformation data",
    long_description=README,
    long_description_content_type="text/markdown",
    scripts = ['src/stripenn.py'],
    author = 'Sora Yoon',
    author_email = 'sora.yoon@pennmedicine.upenn.edu',
    url = 'https://github.com/ysora/stripenn-python',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=False,
    install_requires=["cooler","pandas","numpy","matplotlib","opencv-python","scikit-image","scipy","joblib","tqdm","Bottleneck"],
    entry_points={
    	"console_script": [
    		"stripenn=src.stripenn:main"
    	]
    },
)
