# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("requrements.txt") as f:
    required = f.read().splitlines()

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

"""with open('requirements.txt') as f:
    requirements = f.read()
"""

setup(
    name='celloracle_streamlit',
    version="0.0.1",
    description='celloracle_utility_functions',
    long_description=readme,
    python_requires='>=3.6',
    classifiers=[# How mature is this project? Common values are
                #   3 - Alpha
                #   4 - Beta
                #   5 - Production/Stable
                'Development Status :: 4 - Beta',

                # Indicate who your project is intended for
                'Intended Audience :: Developers',
                'Topic :: Software Development :: Build Tools',

                # Pick your license as you wish (should match "license" above)
                # 'License :: OSI Approved :: MIT License',

                # Specify the Python versions you support here. In particular, ensure
                # that you indicate whether you support Python 2, Python 3 or both.
                'Programming Language :: Python :: 3.6',
                'Programming Language :: Python :: 3.7'
            ],
    install_requires=required,
    author='Kenji Kamimoto at Samantha Morris Lab',
    author_email='kamimoto@wustl.edu',
    #url='https://github.com/morris-lab/CellOracle',
    license=license,
    package_data={"celloracle": []},
    packages=["celloracle_streamlit",
    "celloracle_streamlit.utility",
    "celloracle_streamlit.streamlit_network_analysis_app",
    "celloracle_streamlit.streamlit_perturb_simulation_app",
     "celloracle_streamlit.applications", "celloracle_streamlit.visualizations"],
    #entry_points={'console_scripts':['seuratToAnndata = celloracle.data_conversion.process_seurat_object:main']}

)
