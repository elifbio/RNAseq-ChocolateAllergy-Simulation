#Elif_ChocolateAllergy_RNAseq/

#Transcriptomic Signature of Chocalate Allergy: Simulated RNA-Seq



#Import required packages
import numpy as np
import pandas as pd

#Simulated gene expression counts

# Simülasyon için sabit rastgelelik
np.random.seed(42)

# Gen isimleri (istersen kişisel bağlam ekleyebiliriz)
genes = ['IL4', 'IL13', 'FCER1A', 'HDC', 'CCL2', 'GENEX']


# Kontrol ve alerji koşulları için gen ifadeleri (poisson dağılımı)

# lam değerlerini 6 gen x 3 örnek olacak şekilde yay
lam_control = np.tile([50, 60, 30, 20, 40, 100], (3,1)).T  # shape (6,3)
lam_allergy = np.tile([90, 100, 80, 65, 70, 130], (3,1)).T  # shape (6,3)

# Poisson dağılımı ile sayım verisi üret
control = np.random.poisson(lam=lam_control)
allergy = np.random.poisson(lam=lam_allergy)

# Sayım tablosunu oluştur
counts = pd.DataFrame(
    data = np.hstack([control, allergy]),
    index = genes,
    columns = ['Ctrl_1', 'Ctrl_2', 'Ctrl_3', 'Allergy_1', 'Allergy_2', 'Allergy_3']
)

counts

print(counts)




#Data transformation and statistical analyses




#Visualization



#Biological Interpretation




#Final insights and future description