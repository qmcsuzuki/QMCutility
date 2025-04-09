from setuptools import setup, find_packages

setup(
    name="QMCutility",  # Name of the package
    version="0.1.05",  # Version number
    packages=find_packages(),  # Automatically find packages in the repository
    install_requires=[  # List any dependencies here, if applicable
        "numpy",  # Example dependency
    ],
    # Optional metadata
    author="Kosuke Suzuki",
    author_email="your.email@example.com",
    description="A utility for Quasi Monte Carlo simulations",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/qmcsuzuki/QMCutility",  # Replace with your repository URL
)
