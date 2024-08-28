### Comprend pas pq même avec le bon flux ça marche pas...
    'Cent': 0.6,
    'Cdet': 1.99,       # 'Cdet': 2.5,
    'wp_a': 1.5,
    'wp_b': 1.5,      # 'wp_b': 1.
    'wp_bp': 0.003*250,     #      0.002,
    'up_c': 0.5,
    'vp_c': 0.5,
    'bc_ap': 0.2,    #0.3,
    'delta_bkg': 0.02*250,   # 0.02,
    'output_filename': 'run',
    'wp0':-1.e-08,
    'write_netcdf': True
    + pas de modif de ent/det

![alt text](image.png)
![alt text](image-1.png)

### je pense qu'il doit y avoir échange entre u et v (et wu, wv) dans MNH...
- Ca expliquerait aussi pq les panaches partaient pas dans le bon sens...
- MAIS d'après calcul d'ondes inertielles (cf MAthsFluids): $\underline{U}_m(t)=\frac{u_*^2}{h(t)if} (1-e^{-ift})$ + plot sur wolfram c'est les bons signes pour W05_C500_CORIO...
- [X] vérifier direction sur NTrad drift LES vs litté : cohérent avec Stranéo (=pas d'erreur ds MNH)  
- [ ] demander à Florian de faire tourner WANG1 pour vérifier...
- [W] vérifier signe de la CL cohérente de TKE quand wp0 est fort... --> fait dans scm.class, ne change rien...


![alt text](image-2.png)
![alt text](image-3.png)

    'Cent': 0.3,    
    'Cdet': 1.99,       # 'Cdet': 2.5,
    'wp_a': 1.,
    'wp_b': 1.,      # 'wp_b': 1.
    'wp_bp': 0.003*250,     #      0.002,
    'up_c': 0.5,
    'vp_c': 0.5,
    'bc_ap': 0.35,    #0.3,
    'delta_bkg': 0.02*250,   # 0.02,
    'output_filename': 'run',
    'wp0':-1.e-02,

### Pb de CL de U ? 

![alt text](image-4.png)
cf wu coincident proche surface pour les et scm, mais dz_u ne part pas dans le même sens:
- [ ] changer à la main dans le fortran pour voir
    - dans advance_dyn_mf -uFlx ne fait rien
    - dans advance_dyn_mf changer + FC en - FC fait des instabilités pas bonnes
    - 
- [ ] demander à Florian simu CROCO, pour voir si c'est pas un pb de CL de MNH...

### explication delta0 plus fort
doit ere du au fait que le terme c'est $\delta_0 / h$. Mais avec rot, il faut remplacer $h$ par le rayon $r$ qui doit tendre vers $l_{\rm rot}$ quand $Ro \ll 1$ ? Utiliser une interpolation en tanh entre lrot et h? Besoin d'être justifiée par plusieurs LES... Ou même LES a des instanst différents ? + utiliser longueur intégrale de WANG.

Avec une modulation cff = tanh(Ro^0.37) comme dans Wang de Cent, Cdet (réduit), et de lrot = h x cff pour detla0 et wp_bp, on arrive à retrouver bon résultats sans les coeff (sauf ) !!  