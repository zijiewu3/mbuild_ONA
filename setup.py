from setuptools import setup

setup(
    name="mbuild_ona",
    install_requires="mbuild",
    entry_points={"mbuild.plugins": ["DNA = mbuild_ona.mbuild_ona:DNA","DNA_ds = mbuild_ona.mbuild_ona:DNA_ds"]},
    py_modules=["mbuild_ona"],
)
