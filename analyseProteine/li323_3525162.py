# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 12:51:04 2016

@author: Licia AMICHI 
"""

                                    #*********IMPORTS**************#


#Permet d'accéder aux fonctions mathématiques ici l'importe est utilisé pour pouvoir accéder à la fonction log base 2
import math

#Matploatlib permet de faire des graphes #Pyplot permet de rendre le fonctionnement de matplotlib comme celui de Matlab 
import matplotlib.pyplot as plt

#Itemgetter effectue des appels à travers _getitem_() ce qui facilite les tries 
from operator import itemgetter

                                    #**********CONSTANTES************#


#Les contantes utilisées dans le projets 
M=5643                  #Nombre de ligne dans Dtrain
L=48                    #Nombre de colonne dans Dtrain
A=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'] #Liste des acides aminées
q=21                    #Nombre d'acides aminées
N=114                   #Nombre de colonne dans le fichier test_seq
                        #Pour un teste sur un autre fichier il faut modifier cette valeur et la remplacer par le nombre
                        #de colonne de la séquence considérée


                                    #-------------PARTIE CODE-------------------#


                            #*********************PREMIÈRE PARTIE*************************#


#PREMIÈRE FONCTION:
def Occurrence_poid(fichier):
    #La structure utilisée: une liste de dictionnaire
    #chaque dictionnaire correspond à une position
    #liste=[dictionnairePosition1,dictionnairePosition2... ]
    #un dictionnaire a pour clé une acide aminée et pour valeur son nombre d'occurrence dans la colonne considérée
    #dictionnairePosition1={(AcideAminée: nombre d'occurrence dans la colonne 1),(..)...}
    #dictionnairePositionN={(AcideAminée: nombre d'occurrence dans la colonne N),(..)...}

    liste=list()
    with open(fichier,"r") as f:
        for l in f:
            if l[0]!='>' : 
                for j in range(L):
                    d=dict()
                    if len(liste)<=j:
                        d[l[j]]=1
                        liste.append(d)
                    else:
                        d=liste[j]
                        if  l[j] in d.keys() :
                            d[l[j]]=d[l[j]]+1
                        else:
                            d[l[j]]=1
    f.close()

    #On a finalement une matrice de taille M*L qui associe à chaque acide aminée et sa position ---> son nombre d'occurrence
    
    #Calcul du poid en appliquant la formle (3)
    w = dict()
    for i in range(len(liste)):
        d=liste[i]
        for j in A:
            if j in d.keys():
                w[(j,i)]=float(d[j]+1)/float(M+q)
            else:
                w[(j,i)]=float(1)/float(M+q)
                
    return w



#DEUXIÈME FONCTION:
def entropie_relative(w):
    s=dict()
    w1=dict()
    acide_plus_cons=dict()
    #Calcul de l'entropie relative en appliquant la formule (4)
    for i in range(L):
        s[i]=0.0
        for j in A:
            p=math.log(w[(j,i)],2)
            s[i]=s[i]+w[(j,i)]*p
            #W1 est un dictionnaire qui contient une acide aminée et son poid
            w1[j]=w[(j,i)]      
            #utilisé pour extraire l'acides aminée la plus conservée
        s[i]=s[i]+math.log(q,2)
        w2=w1.items()
        #trie de la liste d'acides aminée selon l'ordre décroissant de leurs poids
        w2.sort(key=itemgetter(1),reverse=True)
        l=w2[:1]                #l'acide aminée la plus conservée pour la position i
        #dictionnaire qui contient pour chaque position l'acide aminée la plus conservée
        acide_plus_cons[i]=(l[0][0]) 

    #Les 3 positions les plus conservées
    w2=s.items()
    w2.sort(key=itemgetter(1),reverse=True)
    l=w2[:3]
    print('La 1ere position plus conservée est :'+str(l[0])+' l''acide aminée la plus conservée :'+str(acide_plus_cons[l[0][0]]))
    print('La 2eme position plus conservée est :'+str(l[1])+' l''acide aminée la plus conservée :'+str(acide_plus_cons[l[1][0]]))
    print('La 3eme position plus conservée est :'+str(l[2])+' l''acide aminée la plus conservée :'+str(acide_plus_cons[l[2][0]]))

    #Graphe:
            #axe des abscisses: les positions (de 0 à L-1)
            #axe des ordonnées: l'entropie relative
    plt.plot(s.keys(),s.values())
    plt.xlabel('Position')
    plt.ylabel('Entropie relative')
    plt.title('Entropie relative en fonction des positions')
    plt.show()

    #Création du fichier entropie_relative.txt en faisant appel à la fonction "creer_fichier_entropie_relative"
    creer_fichier_entropie_relative(s,acide_plus_cons)
    return s

def creer_fichier_entropie_relative(entropie,listeacide):
    fichier=open('entropie_relative.txt','w')
    fichier.write('Position:\tEntropie_relative:\t\tArg_max:\n')
    for i in range(L):
        fichier.write(str(i+1)+'\t\t'+str(entropie[i])+'\t\t\t'+str(listeacide[i])+'\n')
    fichier.close()
    
#TROISIÈME FONCTION   
def generer_modele_null(w):   #Application directe de la formule (8)
    mdn=dict()
    for j in A:
        for i in range(L):
            if j in mdn.keys():
                mdn[j]=mdn[j]+w[(j,i)]
            else:
                mdn[j]=w[(j,i)]
        mdn[j]=float(mdn[j])/float(L)
    #Retourne un dictionnaire qui associe chaque acide aminée à sa fréquence
    return mdn


#QUATRIÈME FONCTION
def log_vraisemblance(w,fich,mdn):  #Application directe de la formule (9)
    with open(fich,"r") as f:
        for l in f:
            if l[0]!='>' :
                lb=dict()
                for k in range(N-L):
                    seq=l[k:k+L]
                    for i in range (L):
                        if k in lb.keys():
                            lb[k]=lb[k]+math.log((w[(seq[i],i)]/mdn[seq[i]]),2)
                        else:
                            lb[k]=math.log((w[(seq[i],i)]/mdn[seq[i]]),2)
                    if lb[k]>0 :
                        #Impression des sequences probable dans le modèle considéré (ici le modèle nul)
                        print('\nLa séquence qui débute à l''indice: '+str(k)+' a pour log vraisemblance: '+str(lb[k])+'\nAppartient à la famille Dtrain\n\n'+str(seq))
    f.close()

    #Création du fichier log_vraisemblance.txt en faisant appel à la fonction "creer_fichier_log_vraisemblance"
    creer_fichier_log_vraisemblance(lb)
    #Graphe
            #axe des abscisses: indice du début de la sous-séquence (de 0 à N-L)
            #axe des ordonnées: la log_vraisemblance correspondante
    plt.plot(lb.keys(),lb.values())
    plt.show()
    return lb

def creer_fichier_log_vraisemblance(lb):
    fichier=open('log_vraisemblance.txt','w')
    fichier.write('Position:\tlog_vraisemblance:\n')
    for i in range(N-L):
        fichier.write(str(i+1)+'\t\t'+str(lb[i])+'\n')
    fichier.close()
        


                                #***********************PARTIE DEUX ******************************#

#PREMIÈRE FONCTION: LA MEME QUE PRÉCEDEMMENT


#DEUXIÈME FONCTION:
def Occurence_poid_double_indice(fichier):
    #initialisation de la matrice du nombre d'occurrences à 0
    n=dict()
    for i in range(L):
        for j in range(L):
            for a1 in A:
                for a2 in A:
                    n[(i,j,a1,a2)]=0

    #Mis à jour de la matrice d'occurrences              
    with open(fichier,"r") as f:
        for l in f:
            if l[0]!='>' :
                for i in range(L-1):
                    for j in range(i+1,L):
                        n[(i,j,l[i],l[j])]+=1
    f.close()
    w=dict()
    #Calcul du poid application directe de la formule (11)
    for i in range(L-1):
        for j in range(i+1,L):
            for a1 in A:
                for a2 in A:
                    w[(i,j,a1,a2)]=float(n[(i,j,a1,a2)]+float(1.0/q))/float(M+q)
    
    return w


#TROISIÈME FONCTION:
def information_mutuelle(w1,w2):  #Quantification des corrélations
    m=dict()
    for i in range(L-1):
        for j in range(i+1,L):
            for a1 in A:
                for a2 in A:
                    if (i,j) in m.keys():
                        if w2[(i,j,a1,a2)]>0:
                            m[(i,j)]=m[(i,j)]+w2[(i,j,a1,a2)]*math.log(w2[(i,j,a1,a2)]/(w1[(a1,i)]*w1[(a2,j)]),2)
                    else:
                        if w2[(i,j,a1,a2)]>0:
                            m[(i,j)]=w2[(i,j,a1,a2)]*math.log(w2[(i,j,a1,a2)]/(w1[(a1,i)]*w1[(a2,j)]),2)
    return m

#Lecture du fichier donnée en paramètre qui "distances.txt" et renvoie un dictionnaire où les positions (i,j) forment la clé
#et la distance la valeur associée à la clée
def generer_dictionnaire_distance(fichier):
    dicDistance=dict()
    with open(fichier,"r") as f:
        for l in f:
            l2=l.split()
            dicDistance[(l2[0],l2[1])]=l2[2]
    return dicDistance


#QUATRIÈME FONCTION:
def Comparer_paires_correlees(m,dicDistance):
    #création du fichier paires_de_positions.txt en faisant appel à la fonction "creer_fichier_paires_positions"
    creer_fichier_paires_positions(m)
    #dicCoo=dict()
    resultat=dict()
    w1=m.items()
    w1.sort(key=itemgetter(1),reverse=True)  #Trie des Mij
    cmp=0
    l=w1[:50]
    #Boucle utilisée pour selectionner au fur et à mesur les des paires entre 1 et 50
    for i in range(50):
        a=dicDistance[(str(l[i][0][0]),str(l[i][0][1]))]
        if float(a) <= 8.0:
            cmp=cmp+1
        resultat[i+1]=float(cmp)/float(i+1)

    #Graphe
            #axe des abscisses: nombre de paires selectionnées
            #axe des ordonnées: fraction du nombre de paires selectionnées avec distance <=8
    plt.plot(resultat.keys(),resultat.values())
    plt.xlabel('Nombre de paire')
    plt.ylabel('Fraction')
    plt.title('Fraction en fonction du nombre de paires')
    plt.show()
        
       
def creer_fichier_paires_positions(m):
    fichier=open('paires_de_positions.txt','w')
    fichier.write('Couple (i,j):\tM[i,j]:\n')
    for i in range(L-1):
        for j in range(i+1,L):
            fichier.write('('+str(i+1)+','+str(j+1)+')\t\t'+str(m[(i,j)])+'\n')
    fichier.close()
        
def accueil():
    print('------------------------------------------------------------------------------------------------------------------------')
    print('\t\t\t\t\t**Statistique en Bioinformatique**\n')
    print('\t\t\t\t  *Analyse statistique d''une famille de pretéines*\n')
    print('------------------------------------------------------------------------------------------------------------------------')
    print('\n')


#FOONCTION MAIN QUI ENGLOBE LES APPELS
def main():
    accueil()
    w1=Occurrence_poid('Dtrain.txt')
    print('I- Entropie relative:\n')
    s=entropie_relative(w1)

                           #Tests donnés
    print('\n-------------------------Tests-----------------------------')
    print('W0(''-'')='+str(w1[('-',0)]))
    print('S0='+str(s[0]))

    print('\nII- La log-vraisemblance:\n')
    l=log_vraisemblance(w1,'test_seq.txt',generer_modele_null(w1))

                           #Tests donnés
    print('\n-------------------------Tests-----------------------------')
    print('l(b0,....,b(L-1))='+str(l[0]))


    w2=Occurence_poid_double_indice('Dtrain.txt')
    m=information_mutuelle(w1,w2)
    dico=generer_dictionnaire_distance('distances.txt')
    print('\n\nIII- Co-evolution:')

                            #Tests donnés
    print('\n-------------------------Tests-----------------------------')
    print('M(0,1)='+str(m[(0,1)]))
    print('\n\n')
    Comparer_paires_correlees(m,dico)
main()



    
