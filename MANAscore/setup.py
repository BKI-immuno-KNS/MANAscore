from setuptools import setup, find_packages

setup(
    name="MANAscore",
    version="1.0.1",
    description="A package for prediction the MANAscore for TIL",
    author="Zhen Zeng",
    author_email="zzeng22@jhmi.edu",
    url="https://github.com/BKI-immuno-KNS/MANAscore/MANAscore",  
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scikit-learn==1.4.2",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            'manascore=MANAscore:main',  # Adjust if 'main' is located elsewhere
        ],
    },
)
