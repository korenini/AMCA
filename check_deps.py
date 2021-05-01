import subprocess


def adc():
    adversion = subprocess.run(["asciidoc", "--version"], capture_output=True)
    
    if adversion.returncode == 0:
        ad_ok = True
    else:
        ad_ok = False
        print("Asciidoc not installed.")
    
    return(ad_ok)


def check_mpl():
    try:
        import matplotlib
        mpl = True
    except ModuleNotFoundError:
        mpl = False
        print("Matplotlib not installed.")

    return mpl
