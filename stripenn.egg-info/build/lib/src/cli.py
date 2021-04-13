import typer
from src import stripenn
from src import getStripe
from src import seeimage
from src import score
import multiprocessing

app = typer.Typer()

@app.command('compute')
def execute(
    cool: str = typer.Option(..., "--cool",help="Path to cool file"),
    out: str = typer.Option(..., "--out", "-o", help="Path to output directory"),
    norm: str = typer.Option('KR',"--norm",help="Normalization method. It should be one of the column name of Cooler.bin(). Check it with Cooler.bins().columns (e.g., KR, VC, VC_SQRT)"),
    chrom: str = typer.Option('all', "--chrom", "-k", help="Set of chromosomes. e.g., 'chr1,chr2,chr3', 'all' will generate stripes from all chromosomes"),
    canny: float = typer.Option(2.5, "--canny", "-c", help="Canny edge detection parameter."),
    minL: int = typer.Option(10,'--minL','-l', help="Minimum length of stripe."),
    maxW: int = typer.Option(8, '--maxW','-w', help="Maximum width of stripe."),
    maxpixel: str = typer.Option('0.95,0.96,0.97,0.98,0.99','--maxpixel','-m', help="Percentiles of the contact frequency data to saturate the image. Separated by comma"),
    numcores: int = typer.Option(multiprocessing.cpu_count(), '--numcores','-n', help='The number of cores will be used.'),
    pvalue: float = typer.Option(0.1,  '--pvalue','-p', help='P-value cutoff for stripe.'),
    mask: str = typer.Option('0', "--mask", help='Column coordinates to be masked. e.g., chr9:12345678-12345789'),
    slow: bool= typer.Option(False,"-s", help='Use if system memory is low.')
):
    """Finds stripe coordinates from 3D genomic data
    """
    stripenn.compute(cool, out, norm, chrom, canny, minL, maxW, maxpixel, numcores, pvalue, mask, slow)

@app.command('seeimage')
def seeimag(
        cool: str = typer.Option(..., "--cool",help="Path to cool file"),
        position: str = typer.Option(..., "--position",'-p', help="Genomic position (e.g., chr1:135010000-136000000)"),
        maxpixel: str = typer.Option('0.95,0.96,0.97,0.98,0.99',"--maxpixel",'-m', help="Quantile for the pixel saturation. (e.g., 0.95)"),
        out: str = typer.Option('./heatmap.png', "--out", "-o", help="Path to output directory"),
        norm: str = typer.Option('KR',"--norm",help="Normalization method. It should be one of the column name of Cooler.bin(). Check it with Cooler.bins().columns (e.g., KR, VC, VC_SQRT)"),
        slow: bool= typer.Option(False,'-s' , help='Use if system memory is low.')
):
    """ Draws heatmap image of given position and color saturation parameter (maxpixel).
    """
    seeimage.seeimage(cool, position, maxpixel, norm, out, slow)
    return 0

@app.command('score')
def scoring(
    cool: str = typer.Option(..., "--cool",help="Path to cool file"),
    coordinates: str = typer.Option(..., "--coord",'-c', help="Path to stripe coordinate table"),
    norm: str = typer.Option('KR',"--norm",help="Normalization method. It should be one of the column name of Cooler.bin(). Check it with Cooler.bins().columns (e.g., KR, VC, VC_SQRT)"),
    numcores: int = typer.Option(multiprocessing.cpu_count(), '-n','--numcores', help='The number of cores will be used.'),
    out: str = typer.Option('scores.out','--out','-o',help='Path to output file'),
    mask: str = typer.Option('0', "--mask", help='Column coordinates to be masked. e.g., chr9:12345678-12345789')
):
    """ Calculates p-value and stripiness of given stripes based on given 3D genome conformation data.
    """
    score.getScore(cool, coordinates,norm,numcores,out,mask)


def main():
    app()

if __name__ == "__main__":
    app()
