import argparse
import matplotlib.pyplot as plt

IMAGE_RESOLUTION = (2000, 2000)  # Résolution de l'image PNG

def plot_haplotype_matrix(file_path, output_file):
	# Lecture du fichier
	with open(file_path, 'r') as file:
		lines = file.readlines()[1:]  # Ignorer la première ligne

	# Création de la matrice à partir des lignes lues et inversion de l'ordre des lignes
	matrix = [[int(char) for char in line.strip()] for line in lines][::-1]

	# Création du plot
	plt.figure(figsize=(IMAGE_RESOLUTION[0]/100, IMAGE_RESOLUTION[1]/100))
	plt.imshow(matrix, cmap='binary', aspect='auto', origin='lower')
	plt.axis('off')  # Supprimer les axes
	plt.xticks([])   # Supprimer les ticks sur l'axe x
	plt.yticks([])   # Supprimer les ticks sur l'axe y
	plt.tight_layout(pad=0)  # Supprimer les marges

	# Sauvegarde de l'image dans un fichier PNG avec la résolution spécifiée
	plt.savefig(output_file, bbox_inches='tight', pad_inches=0, dpi=IMAGE_RESOLUTION[0]/(IMAGE_RESOLUTION[0]/100))
	plt.close()  # Fermer la figure après l'exportation

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Generate a haplotype matrix image from a msms file produced by msms')
	parser.add_argument('--input', type=str, default='simulation.msms', help='path to the input file (default: %(default)s)')
	parser.add_argument('--output', type=str, default='simulation.png', help='path to save the output PNG image (default: %(default)s)')

	args = parser.parse_args()
	plot_haplotype_matrix(args.input, args.output)
