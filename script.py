import sys
import numpy as np    #do macierzy kontaktów
from Bio import PDB     #do analizy struktury białka z pliku PDB.
import matplotlib.pyplot as plt  #by wykres mieć

def save_plot_contact_map(contact_map, output_image):
    plt.imshow(contact_map, cmap='viridis', interpolation='none') #brak interpolacji, kolory virdis
    plt.title('Contact Map')   #tytul
    plt.xlabel('Residue Index')  #osX
    plt.ylabel('Residue Index')  #osY
    plt.colorbar()   #tworzy legendę kolorów, która przyporządkowuje kolorystykę wykresu do odpowiadających wartości w mapie kolorów
    plt.savefig(output_image, format='pdf')  # Zapisujemy wykres do pliku PDF
#
def calculate_contact_map(structure, threshold=8.0):
    # przekształca generator w listę i uzyskaj liczbę reszt
    # residues = list(structure[0].get_residues())
    # num_residues = len(residues)

    # numpy i  macierz kontaktów

    coords = []
    # iteracja przez atomy CA
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] != ' ':
                    continue
                ca = residue['CA']
                coord = ca.get_coord()
                coords.append(coord)



    num_atoms = len(coords)
    contact_map = np.zeros((num_atoms, num_atoms), dtype=int)

    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            distance = np.linalg.norm(coords[i] - coords[j])
            if distance < threshold:
                contact_map[i][j] = 1
                contact_map[j][i] = 1


    return contact_map





def plot_contact_map(contact_map):
    #tu wykorzystamy matplotliba do stworzenia wizualizacji mapy kontaktów za pomocą wykresu typu heatmap
    plt.imshow(contact_map, cmap='viridis', interpolation='none') #funkcji imshow z biblioteki Matplotlib do wyświetlenia macierzy kontaktów jako obrazu,  kolory z palety virdis, brak interpolacji
    plt.title('Contact Map') #tytuł
    plt.xlabel('Residue Index') # oś X odnosi się do indeksów reszt białka w macierzy kontaktów.
    plt.ylabel('Residue Index')  # oś Y  odnosi się do indeksów reszt białka w macierzy kontaktów.
    # plt.colorbar()  # Dodaje pasek kolorów, który przyporządkowuje kolorystykę z wykresu do odpowiadających im wartości w macierzy kontaktów.
    plt.show()   #wyświetla

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]

    # Wczytujemy strukturę białka
    pdb_parser = PDB.PDBParser(QUIET=True)  #tworzymy obiekt parsera plików PDB z modułu Biopython Bio.PDB Argument QUIET=True wskazuje, że parser działa w trybie cichym czyli nie wyświetla ostrzeżeń ani komunikatów odczytu podczas przetwarzania pliku.
    structure = pdb_parser.get_structure('protein', pdb_file)  #wykorzystujemy parser PDB do odczytania struktury białka z pliku PDB.

    # Wyznaczamy mapę kontaktów
    contact_map = calculate_contact_map(structure)

    # Zapisujemy wykres mapy kontaktów do pliku
    save_plot_contact_map(contact_map, 'contactMap.pdf')

    # Wyświetlamy mapę kontaktów
    plot_contact_map(contact_map)
