from setuptools import setup
from setuptools import find_packages


setup(name='fargene',
        version='0.1',
        description='',
        url='',
        author='Fanny Berglund',
        license='MIT',
        packages=find_packages(),
        include_package_data=True,
        install_requires=['matplotlib<=3.5','numpy<=1.21'],
        entry_points={
            'console_scripts': [
                'fargene=fargene_analysis.fargene_analysis:main',
                'pick_long_reads=fargene_analysis.pick_long_reads:main',
                'fargene_model_creation=fargene_model_creation.create_and_optimize_model:main',
                ],
            },
        zip_safe=False)


