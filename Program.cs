using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.Threading;
using System.Drawing;

namespace InitailLife
{
    class Program
    {
        static string[] extraire_parametres(string chemin_param)
        {
            string[] lignes = System.IO.File.ReadAllLines(chemin_param);
            string[] retour = new string[lignes.Length];
            for(int i=0; i<lignes.Length;i++)
            {
                string[] spt = lignes[i].Split('=');
                retour[i] = spt[1];
            }
            return retour;
        }
        static void Main(string[] args)
        {
            string chemin_param = "D:/EvoV4/Parametres.txt";

            string param_prerempli = "Creer_nouvel_univers(bool)=+'\n'+Frame_a_charger(int)=+'\n'+Conserver_contraintes_environementales(bool)=+'\n'+Sauver_simulation(bool)=+'\n'+Creer_rendu(bool)=+'\n'+Chemin_rendu(string)=+'\n'+Chemin_simulation(string)=+'\n'+chemin_log(string)=+'\n'+Nb_de_frames_a_simuler(int)=+'\n'+Dimention_x_image_sortie(int)=+'\n'+Dimention_y_image_sortie(int)=+'\n'+Seed(int)=+'\n'+dt(double)=+'\n'+Taille_d_une_case_en_metres(double)=+'\n'+Nombre_de_cases_par_cote(int)=+'\n'+Epaisseur_Cerveau(int)=+'\n'+Nb_Neurones_Couche_1(int)=+'\n'+Variabilite_nb_neurones_couche_1(int)=+'\n'+...//A_remplacer_par_les_couches_de_neurones_suivantes_puis_supprimer_cette_ligne+'\n'+Max_attaque(double)=+'\n'+Max_horloge_interne(int)=+'\n'+Max_defense(double)=+'\n'+Max_part_don_energie_enfant_01(double)=+'\n'+Max_energie_depart(double)=+'\n'+Max_guerison(double)=+'\n'+Variabilite_regime(double)=+'\n'+Max_temperature_ideale(double)=+'\n'+Variabilite_lumiere_ideale(double)=+'\n'+Max_vitesse(double)=+'\n'+Max_vitesse_rotation(double)=+'\n'+Max_sante_01(double)=+'\n'+Max_altitude(double)=+'\n'+Max_couleur_terrain_nu_01(double)=+'\n'+Max_couleur_vegetation_01(double)=+'\n'+Max_fertilite(double)=+'\n'+Max_sensibilite_vegetation_a_la_temp(double)=+'\n'+Max_temperature_ideale_vegetation(double)=+'\n'+Max_impraticabilite(double)=+'\n'+Max_vegetation(double)=+'\n'+Max_pouvoir_nutritif_vegetation(double)=+'\n'+Max_encombrement_vegetation(double)=+'\n'+Nb_creatures_initiales(double)=+'\n'+Pourriture_influence_temp(double)=+'\n'+Taux_degradation_viande(double)=+'\n'+Croissance_vegetaux(double)=+'\n'+Chauffage_mutuel(double)=+'\n'+Duree_du_jour(double)=+'\n'+Variabilite_vent(double)=+'\n'+Variabilite_temp(double)=+'\n'+Obstruction_relief(double)=+'\n'+Facteur_camouflage(double)=+'\n'+Seuil_mort_vegetaux(double)=+'\n'+Rayon_d_interaction_entre_creatures(double)=+'\n'+Capacite_de_coloration_vegetaux(double)=+'\n'+Gravite(double)=+'\n'+Distance_de_vue(double)=+'\n'+Cout_action_manger(double)=+'\n'+Rendement_reproduction_sexuee(double)=+'\n'+Rendement_clonage(double)=+'\n'+Cout_mouvement(double)=+'\n'+Cout_soin(double)=+'\n'+Cout_temperature(double); Cout_subsistance(double)=+'\n'+Cout_attaque(double)=+'\n'+Masse_par_unite_d_energie(double)=+'\n'+Exposant_calcul_effort(double) Isolation_graisse(double)=+'\n'+Densite_graisse(double)=+'\n'+Densite_attributs(double)=+'\n'+Masse_attributs(double)=+'\n'+Max_viande_mangee_par_kg(double)=+'\n'+Max_vegetation_mangee_par_kg(double)=+'\n'+Frequence_apparition_nouvel_individu(double)=+'\n'+Proba_changement_nom_de_famille_01(double)=+'\n'+Resolution_du_raycast(int)=+'\n'+Taux_mutation_clonage(double)=+'\n'+Taux_mutation_reproduction_sexuee(double)=+'\n'+Frequence_max_combats(double)=+'\n'+Frequence_max_reproduction_sexuee(double)=+'\n'+Frequence_max_clonage(double)=+'\n'+StabiliteOrganesPerception(double)=+'\n'+VariabiliteLuciditeInitiale(double)=+'\n'+Temps(double)=+'\n'+Temperature(double)=";
            
            bool fichier_param_vide = false;
            if (!System.IO.File.Exists(chemin_param))
            {
                fichier_param_vide = true;
                Console.WriteLine("Veuillez remplir les paramètres dans le fichier :"+chemin_param);
            }
            else
            {
                string[] parametres = extraire_parametres(chemin_param);
                //a charger depuis un fichier de config
                System.Runtime.Serialization.IFormatter formatter = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
                System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
                //Definition des parametres
                bool Create = Convert.ToBoolean(parametres[0]);
                int i0 = Convert.ToInt32(parametres[1]);
                bool conserver_environement = Convert.ToBoolean(parametres[2]);
                bool Save = Convert.ToBoolean(parametres[3]);
                bool rendu = Convert.ToBoolean(parametres[4]);
                string Path_rendu =parametres[5] ;
                string Path_simu = parametres[6];
                string chemin_log = parametres[7];
                int NumOfTurns = Convert.ToInt32(parametres[8]);
                int OutputDimX = Convert.ToInt32(parametres[9]);
                int OutputDimY = Convert.ToInt32(parametres[10]);
                int seed = Convert.ToInt32(parametres[11]);
                double Dt = Convert.ToDouble(parametres[12]);
                double Scale = Convert.ToDouble(parametres[13]);//TAILLE UNE CASE
                int TerrainSubdiv = Convert.ToInt32(parametres[14]);
                int EpaisseurCerveaux = Convert.ToInt32(parametres[15]);
                int[,] LayersStructure = new int[EpaisseurCerveaux, 2];
                //Nombre De Neurones              Incertitude
                int n = 16;
                for (int i = 0; i <EpaisseurCerveaux; i ++)
                {
                    LayersStructure[i, 0] = Convert.ToInt32(parametres[n]); LayersStructure[i, 1] = Convert.ToInt32(parametres[n+1]);
                    n += 2;
                }
                Univers U;
                //Valeurs maximales des attributs des creatures
                Creature MaxValues = new Creature();
                MaxValues.Attaque = Convert.ToDouble(parametres[n+0]);
                MaxValues.InternalClock = Convert.ToInt32(parametres[n + 1]);
                MaxValues.Defense = Convert.ToDouble(parametres[n + 2]);
                MaxValues.DonNRJEnfant = Convert.ToDouble(parametres[n + 3]);
                MaxValues.Energie = Convert.ToDouble(parametres[n + 4]);
                MaxValues.Guerison = Convert.ToDouble(parametres[n + 5]);
                MaxValues.Regime = Convert.ToDouble(parametres[n + 6]);
                MaxValues.IdealTemp = Convert.ToDouble(parametres[n + 7]);
                MaxValues.LumiereIdeale = Convert.ToDouble(parametres[n + 8]);
                MaxValues.MaxVitesse = Convert.ToDouble(parametres[n + 9]);
                MaxValues.MaxVitesseRotation = Convert.ToDouble(parametres[n + 10]);
                MaxValues.Sante = Convert.ToDouble(parametres[n + 11]);
                //Valeurs maximales des attributs du terrain
                ElementTerrain Max = new ElementTerrain();
                Max.Altitude = Convert.ToDouble(parametres[n + 12]);
                Max.CouleurTerrainNu = Convert.ToDouble(parametres[n + 13]);
                Max.CouleurVegetation = Convert.ToDouble(parametres[n + 14]);
                Max.Fertilité = Convert.ToDouble(parametres[n + 15]);
                Max.VegetationTempSensibilite = Convert.ToDouble(parametres[n + 16]);
                Max.IdealTempVegetation = Convert.ToDouble(parametres[n + 17]);
                Max.Impraticabilité = Convert.ToDouble(parametres[n + 18]);
                Max.Vegetation = Convert.ToDouble(parametres[n + 19]);
                Max.VegetationPouvoirNutritif = Convert.ToDouble(parametres[n + 20]);
                Max.VegetationEncombrement = Convert.ToDouble(parametres[n + 21]);
                //Constantes environementales
                ConstantesEnvironementales CE = new ConstantesEnvironementales(Convert.ToInt32(parametres[n + 22]), Convert.ToDouble(parametres[n + 23]), Convert.ToDouble(parametres[n + 24]), Convert.ToDouble(parametres[n + 25]), Convert.ToDouble(parametres[n + 26]), Convert.ToDouble(parametres[n + 27]), Convert.ToDouble(parametres[n + 28]), Convert.ToDouble(parametres[n + 29]), Convert.ToDouble(parametres[n + 30]), Convert.ToDouble(parametres[n + 31]), Convert.ToDouble(parametres[n + 32]), Convert.ToDouble(parametres[n + 33]), Convert.ToDouble(parametres[n + 34]), Convert.ToDouble(parametres[n + 35]), Convert.ToDouble(parametres[n + 36]));
                //Constantes Biologiques
                ConstantesBiologiques CB = new ConstantesBiologiques(Convert.ToDouble(parametres[n + 37]), Convert.ToDouble(parametres[n + 38]), Convert.ToDouble(parametres[n + 39]), Convert.ToDouble(parametres[n + 40]), Convert.ToDouble(parametres[n + 41]), Convert.ToDouble(parametres[n + 42]), Convert.ToDouble(parametres[n + 43]), Convert.ToDouble(parametres[n + 44]), Convert.ToDouble(parametres[n + 45]), Convert.ToDouble(parametres[n + 46]), Convert.ToDouble(parametres[n + 47]), Convert.ToDouble(parametres[n + 48]), Convert.ToDouble(parametres[n + 49]), Convert.ToDouble(parametres[n + 50]), Convert.ToDouble(parametres[n + 51]), Convert.ToDouble(parametres[n + 52]), Convert.ToDouble(parametres[n + 61]));
                //Gestion probabilites
                GestionProbabilites GP = new GestionProbabilites(Convert.ToDouble(parametres[n + 53]), Convert.ToDouble(parametres[n + 54]), Convert.ToInt32(parametres[n + 55]), Convert.ToDouble(parametres[n + 56]), Convert.ToDouble(parametres[n + 57]), Convert.ToDouble(parametres[n + 58]), Convert.ToDouble(parametres[n + 59]), Convert.ToDouble(parametres[n + 60]),Convert.ToDouble(parametres[n + 62]));
                //Conditions initiales
                VariablesEnvironementales VE = new VariablesEnvironementales(Convert.ToDouble(parametres[n + 63]), 0, Convert.ToDouble(parametres[n + 64]), 0);
                //Creation d'un nouvel univers et rechargement
                string LOG = "";
                if (Create)
                {
                    if (System.IO.File.Exists(chemin_log))
                    {
                        System.IO.File.Delete(chemin_log);
                    }
                    if (!System.IO.File.Exists(chemin_log))
                    {
                        System.IO.File.Create(chemin_log);
                    }
                    i0 = -1;
                    sw.Reset();
                    sw.Start();

                    U = new Univers(Dt, seed, Scale, TerrainSubdiv, Max, CE, VE, MaxValues, CB, GP, LayersStructure);
                    sw.Stop();
                    LOG += ("--------->Créer univers " + sw.Elapsed + "\n");
                }
                else
                {
                    if (!System.IO.File.Exists(chemin_log))
                    {
                        System.IO.File.Create(chemin_log);
                    }
                    sw.Reset();
                    sw.Start();
                    System.IO.Stream stream = new System.IO.FileStream(Path_simu + i0 + ".bin", System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read);
                    U = (Univers)formatter.Deserialize(stream);
                    stream.Close();
                    if (!conserver_environement)
                    {
                        U.CE = CE;
                        U.GP = GP;
                        U.VE = VE;
                        U.CB = CB;
                    }
                    sw.Stop();
                    LOG += ("--------->Charger Univers " + sw.Elapsed + "\n");
                }
                //Simulation et rendu
                Bitmap b = new Bitmap(OutputDimX, OutputDimY);
                for (int i = i0 + 1; i < NumOfTurns; i++)
                {
                    sw.Reset();
                    sw.Start();
                    U.TourSuivant(MaxValues);
                    LOG += U.log;
                    sw.Stop();
                    sw.Reset();
                    sw.Start();
                    if (rendu)
                    {
                        U.Rendu(OutputDimX, OutputDimY, true, ref b);
                        b.Save(Path_rendu + "_frame_" + i + ".bmp");
                    }
                    if (Save)
                    {
                        System.IO.Stream stream = new System.IO.FileStream(Path_simu + i + ".bin", System.IO.FileMode.Create, System.IO.FileAccess.Write, System.IO.FileShare.None);
                        formatter.Serialize(stream, U);
                        stream.Close();
                    }
                    sw.Stop();
                    if (rendu)
                    {
                        LOG += ("---------" + "Rendu " + sw.Elapsed + "\n");
                    }
                    sw.Reset();
                    LOG += (">---------------------<" + "\n");
                    Console.WriteLine(LOG);
                    using (System.IO.StreamWriter sr = System.IO.File.AppendText(chemin_log))
                    {
                        sr.WriteLine(LOG);
                    }
                    LOG = "";
                }

            }
            if(fichier_param_vide)
            {
                    System.IO.File.WriteAllText(chemin_param,param_prerempli);
            }
        }
    }
    [Serializable]public struct Creature
    {
        //Variables
        public double Sante;
        public double Energie;
        public Complex Vitesse;
        public double VitesseRotation;
        public Complex Position;
        public double Attitude;
        public double VolonteReproductive;
        public double Agresivite;
        public double FoccusLineaire;
        public double FoccusAngulaire;
        
        //Variables determinees
        public double TempsDepuisDerniereReproduction;
        public double TempsDepuisDernierClonage;
        public double Age;
        public double Masse;
        public double TemperaturePercue;
        public double VentPercu;
        public double Effort;
        public double Rayon;
        public Complex Alea;

        //Stats
        public double Defense;
        public double Attaque;
        public double MaxVitesse;
        public double MaxVitesseRotation;
        public double Couleur;
        public double Regime;
        public double Guerison;
        public double IdealTemp;
        public double DonNRJEnfant;
        public double LumiereIdeale;

        //Calculs seulement
        public double MasseGraisse;
        public double MasseAttributs;
        public double RayonGraisse;
        public double RayonAttributs;
        public double BilanEnergetique;


        //History
        public int ID;
        public int P1Id;
        public int P2Id;
        public Nom nom;
        public bool ReproducedThisTurn;
        public int clock;
        public int InternalClock;
        public bool AMange;
        public bool SEstBattu;

        //Neurones
        public double[] Perception; //Indique quels organes percepteurs sont effectifs (0 inefficace, 1 efficace-- arrondi)
        public Cerveau Brain;
        public Complex[] Outputs;
        //0 Vitesse
        //1 Attitude
        //2 VolonteReprodutive
        //3 Agresivite
        //4 Focus Angulaire
        //5 Focus Lineaire
        //6 Manger
        public Complex[] Inputs;
        //Cf createInputs

    }
    [Serializable]public struct ElementTerrain
    {
        public double Fertilité;
        public double Vegetation;
        public double IdealTempVegetation;
        public double VegetationTempSensibilite;
        public double VegetationPouvoirNutritif;
        public double VegetationEncombrement;
        public double Viande;
        public double Impraticabilité;
        public double Altitude;
        public double CouleurTerrainNu;
        public double CouleurVegetation;
        //VariablesApparentes
        public double CapaciteBloquage;
        public double Coloration;
    }
    [Serializable]public class Univers
    {
        public string log = "";

        double PerteThermiques = 0;
        double PertesMouvement = 0;
        double PertesSubsistance = 0;
        double PertesSante = 0;
        double PertesReproduction = 0;
        double PertesCombat = 0;
        double GainAlimentation = 0;
        //Physique
        public double dt;
        //Population
        Creature[] Creatures;
        List<Creature> Enfants = new List<Creature>();
        bool[] EnVie;
        //Environement
        public double Scale; // Longueur en metres d une case
        public int TerrainSubdiv;
        ElementTerrain[,] Terrain;
        ElementTerrain Max_CentreTerrainValues;
        public ConstantesEnvironementales CE;
        public VariablesEnvironementales VE;
        //Biologie
        Creature Max_CentreProprietesCreatures;
        public ConstantesBiologiques CB;
        //Proba
        public GestionProbabilites GP;
        //Network 
        public int[,] CouchesDuReseauEtVariabilite;
        

        public Univers(double dt,int Seed, double scale, int terrainSubdiv, ElementTerrain max_centreTerrainValues, ConstantesEnvironementales cE, VariablesEnvironementales vE, Creature max_centreProprietesCreatures, ConstantesBiologiques cB, GestionProbabilites gP, int[,] couchesDuReseauEtVariabilite)
        {
            this.dt = dt;
            Scale = scale;
            TerrainSubdiv = terrainSubdiv;
            Max_CentreTerrainValues = max_centreTerrainValues;
            CE = cE;
            VE = vE;
            Max_CentreProprietesCreatures = max_centreProprietesCreatures;
            CB = cB;
            GP = gP;
            CouchesDuReseauEtVariabilite = couchesDuReseauEtVariabilite;
            //PremiereGeneration
            GP.RANDOM = new Random(Seed);
            GenerateTerrain(terrainSubdiv, terrainSubdiv, ref GP.RANDOM, Max_CentreTerrainValues);
            GenererPremiersPeuples(CE.PremiersPeuplesNombre, ref GP.RANDOM, Max_CentreProprietesCreatures);
        }

        //Maths et Outils
        double GaussianRandom(ref Random r)
        {
            double r1 = r.NextDouble();
            double r2 = r.NextDouble();
            if(r1==0)
            {
                r1 = 1;
            }
            double k = -Math.Log(r1);
            double rad = 2.0 * Math.PI*r2;
            return 1 + Math.Cos(rad) * k;
        }
        double GaussianBorneRandom(ref Random r)
        {
            double r1 = r.NextDouble();
            if (r1 == 0)
            {
                r1 = 1;
            }
            return Math.Exp(-Math.Pow(Math.Log(-1 + 1.0 / r1), 2));
        }
        double Sig(double x)
        {
            return (1.0 / (1.0 + Math.Exp(-1.0 * x)));
        }
        double SigInverse(double x)
        {
            if (x <= 0 || x >= 1)
            {
                return 0;
            }
            else
            {
                return -Math.Log(-1.0+1.0/x);
            }
        }
        double Normale(double x)
        {
            return Math.Exp(-Math.Pow(x, 2));
        }
        double NormaleInverse(double x)
        {
            if(x<=0)
            {
                return 0;
            }
            else
            {
                return Math.Sqrt(-Math.Log(x));
            }

        }
        double ProbaPoisson(double Frequence)
        {
            //Proba que l'evenement avec une frequence moyenne arrive durant l intervalle dt au moins une fois
            return 1.0 - Math.Exp(-Frequence * dt);
        }

        Complex GetFoccus(int Percepteur, Complex Position)
        {
            
            Complex RelativePosition = Position - Creatures[Percepteur].Position;
            //La distance est exprimee en metres
            double RelativeDist = Distance(Creatures[Percepteur].Position, Position);//en m
            double RelativeAngle = Creatures[Percepteur].Vitesse.Phase-(RelativePosition).Phase;
            double CoeffEffacement = Creatures[Percepteur].FoccusAngulaire * Math.Abs(RelativeAngle) + Creatures[Percepteur].FoccusLineaire * RelativeDist;
            Complex Output = Complex.FromPolarCoordinates(1.0/(1.0 + (CoeffEffacement)), RelativeAngle);
            return Output;
            
        }
        double Obstruction(int Percepteur,Complex Objectif)
        {
            double obs = 0;
            Complex RelativePosition = Objectif - Creatures[Percepteur].Position;
            int Steps = (int)Math.Max(0,Math.Round(RelativePosition.Magnitude * TerrainSubdiv*GP.RaycastResolution*Math.Sqrt(2)));
            double AltitudeDepart = GetTile(Percepteur).Altitude;
            double AltitudeArrivee = 0;
            double Delta = AltitudeArrivee - AltitudeDepart;
            if (Steps != 0)
            {
                for (int i = 1; i < Steps; i++)
                {
                    Complex Ray = ((float)i / (float)Steps) * RelativePosition;
                    Complex Hit = Creatures[Percepteur].Position + Ray;
                    ElementTerrain T = GetTile(Hit);
                    //Altitude de la ligne de mire en ce point
                    double AltitudeMire = AltitudeDepart + ((double)i / Steps) * Delta;
                    //Difference entre l altitude de la ligne de mire et l'altitude du terrain
                    double DeltaAltitude = AltitudeMire - T.Altitude;
                    if (DeltaAltitude >= 0)
                    {
                        //On somme l'obstruction acumulée, proportionelle a la densite d'obstruction rencontree et a la longueur du parcours,(d'ou le facteur obstruction*Scale/scale^2)
                        obs += Math.Exp(-DeltaAltitude * CE.EncombrementRelief) * (T.CapaciteBloquage / Scale);
                    }
                    else
                    {
                        //La ligne de mire est obstruée
                        return 0;
                    }
                }
            }
            else
            {
                return 1;
            }
            //Obstruction moyenne * distance
            obs = (obs/ Steps)*(RelativePosition.Magnitude*TerrainSubdiv*Scale);
            double f = 1.0 / (1.0 + obs);
            if(double.IsNaN(f))
            {
                f = 0;
            }
            return 1.0-f;
        }
        double Camouflage(int Percepteur)
        {
            ElementTerrain T = GetTile(Percepteur);
            //Le camouflage reduit avec le rayon, et la densité d'obstruction (d'ou le facteur obstruction/scale^2)
            double Cam = (1.0/(Scale*Scale))*(T.CapaciteBloquage+CE.FacteurCamouflage*(1-Math.Abs(Creatures[Percepteur].Couleur-T.Coloration)))/Creatures[Percepteur].Rayon;
            double f = 1.0 / (1.0 + Cam);
            return 1.0-f;
        }
        double AccuiteVisuelle(int Percepteur)
        {
            double FacteurLumiere = 1.0 - Math.Abs(Creatures[Percepteur].LumiereIdeale - VE.Lumiere)/2.0;
            return FacteurLumiere;
        }
        double Distance(Complex PositionA, Complex PositionB)
        {
            return (PositionA - PositionB).Magnitude * TerrainSubdiv * Scale;//En metres
        }

        ElementTerrain GetTile(Complex Position)
        {
            Complex K = GetXYofTileAt(Position);
            return Terrain[(int)K.Real,(int)K.Imaginary ];
        }
        Complex GetXYofTileAt(Complex Position)
        {
            int x = (int)Math.Floor(Math.Min((Position).Real, 1.0 - 10.0 * float.Epsilon) * (double)TerrainSubdiv);
            int y = (int)Math.Floor(Math.Min((Position).Imaginary, 1.0 - 10.0 * float.Epsilon) * (double)TerrainSubdiv);
            x = Math.Min(x, TerrainSubdiv - 1);
            y = Math.Min(y, TerrainSubdiv - 1);
            return new Complex(x, y);
        }
        Complex GetPositionOfTile(int x, int y)
        {
            return new Complex((((double)x + 0.5) / (double)Terrain.GetLength(0)), (((double)y + 0.5) / (double)Terrain.GetLength(1)));
        }
        Complex TassementEntree(Complex x)
        {
            Complex k=Complex.FromPolarCoordinates(Math.Min(Math.Pow(x.Magnitude,0.9),1),x.Phase);
            if (double.IsNaN(k.Real)|| double.IsNaN(k.Imaginary) || double.IsInfinity(k.Real) || double.IsInfinity(k.Imaginary))
            {
                ;
                k = 0;
            }
                return k;
        }
        //Depenses
        void UpdateEffort(int Percepteur)
        {
            // L'effort dépend de la vitesse de déplacement souhaitée par la creature, et est proportionel a la masse et à la pente
            //Intesite comprise entre 0 et (masse*Vitessemax)+1;
            ElementTerrain T = GetTile(Creatures[Percepteur].Position);
            double Maintient = (Creatures[Percepteur].Vitesse.Magnitude+Math.Abs(Creatures[Percepteur].VitesseRotation))* T.CapaciteBloquage;
            double Intensite = Math.Pow(Maintient, CB.ExposantEffort)/1000.0;
            //Calcul de la pente
            ElementTerrain Tx = GetTile(Creatures[Percepteur].Position+ new Complex(1.0/(double)TerrainSubdiv,0));
            ElementTerrain Ty = GetTile(Creatures[Percepteur].Position+ new Complex(0,1.0/(double)TerrainSubdiv));
            double dx = (Tx.Altitude - T.Altitude) / Scale;
            double dy = (Ty.Altitude - T.Altitude) / Scale;
            double effm = (dx * Creatures[Percepteur].Vitesse.Real + dy * Creatures[Percepteur].Vitesse.Imaginary)*CE.Gravite;
            //Effort compris entre 0 et +inf
            Creatures[Percepteur].Effort =  Math.Max(Intensite+effm,0) * Creatures[Percepteur].Masse;
            if (double.IsNaN(Creatures[Percepteur].Effort) || double.IsInfinity(Creatures[Percepteur].Effort))
            {
                Creatures[Percepteur].Effort = 0;
                Creatures[Percepteur].Outputs[0] = 0;
            }
        }
        void TemperatureDepenses(int Percepteur)
        {
            double DeltaTemp = Creatures[Percepteur].TemperaturePercue - Creatures[Percepteur].IdealTemp;
            double ResistanceThermique = CB.IsolationGraisse/(4*Math.PI)* ((1.0/ Creatures[Percepteur].RayonAttributs) -(1.0 /(Creatures[Percepteur].Rayon)));
            double DeltaEffectif = 0;
            if(DeltaTemp<0)
            {
                //La creature a froid
                DeltaEffectif = Math.Abs(DeltaTemp) * (1.0 + Creatures[Percepteur].VentPercu);
            }
            else
            {
                //La creature a chaud
                DeltaEffectif = Math.Abs(DeltaTemp)/ (1.0 + Creatures[Percepteur].VentPercu);
            }
            //La perte est proportionelle au delta percu, ainqi qu'a la ResistanceThermique de la creature
            double perte = Math.Max(CB.CoutTemperature*DeltaEffectif/ResistanceThermique,0);
            Creatures[Percepteur].BilanEnergetique -= perte;
            PerteThermiques += perte;
            if (double.IsInfinity(Creatures[Percepteur].BilanEnergetique) || double.IsNaN(perte))
            {
                Creatures[Percepteur].BilanEnergetique = 0;
            }
        }
        void MouvementDepense(int Percepteur)
        {
            double depense = Creatures[Percepteur].Effort * CB.CoutMouvement;
            Creatures[Percepteur].BilanEnergetique -= depense;
            PertesMouvement += depense;
            if(double.IsInfinity(Creatures[Percepteur].BilanEnergetique))
            {
                Creatures[Percepteur].BilanEnergetique = 0;
            }
        }
        void SanteUpdateEtCout(int Percepteur)
        {
            double gainSante = dt*(1.0 - Creatures[Percepteur].Sante) * Sig(Creatures[Percepteur].Guerison) / (1.0 + Creatures[Percepteur].Effort);
            Creatures[Percepteur].Sante += gainSante;
            double j = CB.CoutSoin * gainSante * Creatures[Percepteur].Masse;
            Creatures[Percepteur].BilanEnergetique -= j;
            PertesSante += j;
        }
        void SubsistanceDepense(int Percepteur)
        {
            double j = Creatures[Percepteur].Masse * CB.CoutSubsistance;
            Creatures[Percepteur].BilanEnergetique -= j;
            PertesSubsistance += j;
        }
        void AppliquerBilanEnergetique(int Percepteur)
        {
            Creatures[Percepteur].Energie += Creatures[Percepteur].BilanEnergetique * dt;
            Creatures[Percepteur].BilanEnergetique = 0;
        }
        //Update Etat et variables internes
        void UpdateHorloges(int Percepteur)
        {
            //L'horloge interne prend un tick
            Creatures[Percepteur].clock = (Creatures[Percepteur].clock + 1) % Creatures[Percepteur].InternalClock;
            //
            Creatures[Percepteur].ReproducedThisTurn = false;
            Creatures[Percepteur].Age+=dt;
            Creatures[Percepteur].TempsDepuisDerniereReproduction+= dt;
            Creatures[Percepteur].TempsDepuisDernierClonage += dt;
        }
        void UpdateMasse(int Percepteur)
        {
            //MAsse des attributs de combat + masse de la graisse
            Creatures[Percepteur].MasseGraisse = Creatures[Percepteur].Energie * CB.MasseNRJ;
            Creatures[Percepteur].MasseAttributs = (Creatures[Percepteur].Defense + Creatures[Percepteur].Attaque) * CB.MasseAttributsDefenseOuAttaque;
            Creatures[Percepteur].Masse = Creatures[Percepteur].MasseGraisse + Creatures[Percepteur].MasseAttributs;
        }
        void UpdateRayon(int Percepteur)
        {
            //En m^3
            double VolumeAttributs = Creatures[Percepteur].MasseAttributs / CB.DensiteAttributs;
            double VolumeGraisse = Creatures[Percepteur].MasseGraisse / CB.DensiteGraisse;
            double RayonAttributs = Math.Pow(3.0*VolumeAttributs/(4.0*Math.PI), 1.0 / 3.0);
            double Rayon = Math.Pow(3.0 * (VolumeAttributs + VolumeGraisse) / (4.0 * Math.PI), 1.0 / 3.0);
            //On Considere une boule d'attributs dans une boule de graisse
            //Vgraisse=4/3*pi*(Rattributs+Rgraisse)^3 - Vgraisse
            double RayonGraisse = Rayon-RayonAttributs;
            Creatures[Percepteur].RayonGraisse = RayonGraisse;
            Creatures[Percepteur].RayonAttributs = RayonAttributs;
            //Racine cubique de la masse/densite graisse
            Creatures[Percepteur].Rayon = Rayon;
        }
        void UpdateVentPercu(int Percepteur)
        {
            double vp = 0;
            Complex pos = Creatures[Percepteur].Position;
            if (VE.Vent.Magnitude == 0)
            {
                vp = 0;
            }
            else
            {
                //On trouve l'origine du vent
                Complex OrigineVent = 0;
                double tx = 0;
                double ty = 0;
                //x0+x*t=1
                tx = (1.0 - pos.Real) / VE.Vent.Real;
                //x0+x*t=0
                if (tx<0)
                {
                    tx = (-pos.Real) / VE.Vent.Real;
                }
                //y0+y*t=1
                ty = (1.0 - pos.Imaginary) / VE.Vent.Imaginary;
                //y0+y*t=0
                if (ty < 0)
                {
                    ty = (-pos.Imaginary) / VE.Vent.Imaginary;
                }
                OrigineVent = pos + Math.Min(tx, ty) * VE.Vent;
                //Le relief obstrue le vent jusqu'a la creature
                double obs = Obstruction(Percepteur, OrigineVent);
                if(obs<0)
                {
                    obs = 0;
                }
                 vp= (1.0-obs) * VE.Vent.Magnitude;
            }
            if(double.IsNaN(vp))
            {
                vp = 0;
            }
            Creatures[Percepteur].VentPercu = vp;
        }
        void UpdateTemperaturePercue(int Percepteur)
        {
            //Temperature apportee par les autres individus et par la graisse
            double tp = VE.Temperature;
            //Creatures proches
            for (int i = 0; i < Creatures.GetLongLength(0); i++)
            {
                if (EnVie[i] && i != Percepteur)
                {
                    double RelativeDist = Distance(Creatures[Percepteur].Position, Creatures[i].Position);//en m
                    //La chaleur apportee depend des ratios de la masse et du rayon de la creature voisine avec ceux du percepteur, et est inversement proportionelle a la distance a la puissance Vent
                    double facteurDistance = Math.Pow(RelativeDist + 1.0, -(Creatures[Percepteur].VentPercu + 1));
                    double facteurProtection = (Creatures[i].Rayon/ Creatures[Percepteur].Rayon) * (Creatures[i].IdealTemp/ Creatures[Percepteur].IdealTemp);
                    double tpAdd = facteurDistance*facteurProtection;
                   if( double.IsNaN(tpAdd) || double.IsInfinity(tpAdd))
                    {
                        tpAdd = 0;
                    }
                    tp += tpAdd*CE.ChauffageMutuel;
                }
            }
            Creatures[Percepteur].TemperaturePercue = tp;
        }
        void UpdateVitesse(int Percepteur)
        {
            //Désaxage et intensite du mouvement
            Creatures[Percepteur].VitesseRotation = (Creatures[Percepteur].Outputs[0].Phase / Math.PI) * Creatures[Percepteur].MaxVitesseRotation;
            Complex Rotation = Complex.FromPolarCoordinates(Creatures[Percepteur].Outputs[0].Magnitude, Creatures[Percepteur].Vitesse.Phase+ Creatures[Percepteur].VitesseRotation * dt)  ;
            if (double.IsNaN(Rotation.Real)  || double.IsNaN(Rotation.Imaginary))
            {
                Creatures[Percepteur].Outputs[0] = 1.0;
                Rotation = 0;
            }
            double CB = (1.0 + GetTile(Percepteur).CapaciteBloquage);
            //Vitesse en m/s, depend de la sante et de l encombrement
            Creatures[Percepteur].Vitesse =  Rotation* Creatures[Percepteur].Sante * Creatures[Percepteur].MaxVitesse /CB;
        }
        void UpdateVivant(int Percepteur, bool AffCause)
        {
            if (EnVie[Percepteur])
            {
                if ((Creatures[Percepteur].Sante <= 0 || Creatures[Percepteur].Energie <= 0))
                {
                    //Si la creature n'a plus de sante ou plus d'energie, elle meurt
                    EnVie[Percepteur] = false;
                    int x = (int)Math.Floor((Creatures[Percepteur].Position).Real * (double)TerrainSubdiv);
                    int y = (int)Math.Floor((Creatures[Percepteur].Position).Imaginary * (double)TerrainSubdiv);
                    x = Math.Max(Math.Min(x, TerrainSubdiv - 1),0);
                    y = Math.Max(Math.Min(y, TerrainSubdiv - 1),0);
                    //On ajoute sa viande
                    Terrain[x, y].Viande += Math.Max(Creatures[Percepteur].Energie, 0);
                    Nom nom = Creatures[Percepteur].nom;
                    double Age = Creatures[Percepteur].Age;
                    if (AffCause)
                    {
                        if ((Creatures[Percepteur].Sante <= 0))
                        {
                            log+=("+" + nom.Prenom + " " + nom.NomDeFamille + " est mort de ses blessures à l'age de " + Age+"\n");
                        }
                        else
                        {
                            log+=("+" + nom.Prenom + " " + nom.NomDeFamille + " est mort de faim à l'age de " + Age+"\n");
                        }
                    }
                }
            }
        }
        void UpdateAttitude(int Percepteur)
        {
            //Attitude affichée aux autres individus, entre 0 et 1
            Creatures[Percepteur].Attitude = (Creatures[Percepteur].Outputs[1].Real+1.0)/2.0;
        }
        void UpdateVolonteReproductive(int Percepteur)
        {
            //A quel point l'individu veut se reproduire, comprise entre -inf et + inf
            Creatures[Percepteur].VolonteReproductive = -Math.Log(-1.0 + 2.0 / (Creatures[Percepteur].Outputs[2].Real + 1.0));
        }
        void UpdateAgressivite(int Percepteur)
        {
            //Attitude de l'individu face au combat, entre 0 et 1 (1= tres agressif)
            Creatures[Percepteur].Agresivite = (Creatures[Percepteur].Outputs[3].Real+1.0)/2.0;
        }
        void UpdateVision(int Percepteur)
        {
            //Donne les valeurs des foccus angulaiures et lineaires, entre 0 et +inf
            Creatures[Percepteur].FoccusLineaire = 2.0/(Creatures[Percepteur].Outputs[5].Real+1)-1;
            Creatures[Percepteur].FoccusAngulaire = 2.0 / (Creatures[Percepteur].Outputs[4].Real + 1) - 1;
            if(double.IsInfinity(Creatures[Percepteur].FoccusLineaire)|| double.IsInfinity(Creatures[Percepteur].FoccusAngulaire))
            {
                Creatures[Percepteur].FoccusLineaire = 0;
                Creatures[Percepteur].FoccusAngulaire = 0;
            }

        }
        void UpdateAlimentation(int Percepteur)
        {
            if (Creatures[Percepteur].Outputs[5].Real > 0)
            {
                //La creature mange
                Creatures[Percepteur].AMange = true;
                int x = (int)Math.Floor((Creatures[Percepteur].Position).Real * (double)TerrainSubdiv);
                int y = (int)Math.Floor((Creatures[Percepteur].Position).Imaginary * (double)TerrainSubdiv);
                x = Math.Min(x, TerrainSubdiv - 1);
                y = Math.Min(y, TerrainSubdiv - 1);
                //A quel point l'individu se nourrit
                double ViandeMangee = Math.Min(Terrain[x, y].Viande,  Creatures[Percepteur].Outputs[5].Real*CB.MaxViandeMangeeParKg* Creatures[Percepteur].Masse*dt);
                double HerbeMangee = Math.Min(Terrain[x, y].Vegetation, Creatures[Percepteur].Outputs[5].Real * CB.MaxVegetationMangeeParKg * Creatures[Percepteur].Masse*dt);
                //La nourriture est consommee
                Terrain[x, y].Viande -= ViandeMangee;
                Terrain[x, y].Vegetation -= HerbeMangee;
                //Fournit de l'energie en fonction du regime, et coute de l'energie
                double Bilan = (Creatures[Percepteur].Regime * ViandeMangee + (1.0 - Creatures[Percepteur].Regime) * HerbeMangee* Terrain[x, y].VegetationPouvoirNutritif) - CB.CoutManger * Creatures[Percepteur].Outputs[5].Real*dt;
                //L'energie est gagnee
                double j=Bilan/dt;
                Creatures[Percepteur].BilanEnergetique += j;
                GainAlimentation += j;

            }
            else
            {
                //La creature ne mange pas
                Creatures[Percepteur].AMange = false;
            }
        }
        void CreateInternalInputs(int Percepteur)
        {

            //FeedBack les outputs
            for (int i = 0; i < 7; i++)
            {
                Creatures[Percepteur].Inputs[i] = TassementEntree(Creatures[Percepteur].Outputs[i]);
            }
            //Input les variables environementales + le terrain local
            ElementTerrain T = GetTile(Percepteur);
            Complex LocalHerb = T.Vegetation;
            Complex LocalViande = T.Viande;
            Complex LocalEncombrement = T.CapaciteBloquage;
            Complex LocalAltitude = T.Altitude;
            Complex LocalColoration = T.Coloration;
            Creatures[Percepteur].Inputs[7] = TassementEntree(LocalHerb);
            Creatures[Percepteur].Inputs[8] = TassementEntree(LocalViande);
            Creatures[Percepteur].Inputs[9] = TassementEntree(LocalEncombrement);
            Creatures[Percepteur].Inputs[10] = TassementEntree(VE.Temperature - Creatures[Percepteur].IdealTemp);
            Creatures[Percepteur].Inputs[11] = TassementEntree(VE.Lumiere);
            Creatures[Percepteur].Inputs[12] = TassementEntree(Creatures[Percepteur].VentPercu);
            Creatures[Percepteur].Inputs[51] = TassementEntree(LocalAltitude);
            Creatures[Percepteur].Inputs[53] = TassementEntree(LocalColoration);

            //Input en absolue les variables + le regime
            Creatures[Percepteur].Inputs[13] = TassementEntree(Creatures[Percepteur].Agresivite);
            Creatures[Percepteur].Inputs[14] = TassementEntree(Creatures[Percepteur].Attitude);
            Creatures[Percepteur].Inputs[15] = TassementEntree(Creatures[Percepteur].Age);
            Creatures[Percepteur].Inputs[16] = TassementEntree(Creatures[Percepteur].Effort);
            Creatures[Percepteur].Inputs[17] = TassementEntree(Creatures[Percepteur].FoccusAngulaire);
            Creatures[Percepteur].Inputs[18] = TassementEntree(Creatures[Percepteur].FoccusLineaire);
            Creatures[Percepteur].Inputs[19] = TassementEntree(Creatures[Percepteur].Energie);
            Creatures[Percepteur].Inputs[20] = TassementEntree(Creatures[Percepteur].Masse);
            Creatures[Percepteur].Inputs[21] = TassementEntree(Creatures[Percepteur].Sante);
            Creatures[Percepteur].Inputs[22] = TassementEntree(Creatures[Percepteur].TemperaturePercue);
            //Horloge interne
            Creatures[Percepteur].Inputs[48] = TassementEntree((double)Creatures[Percepteur].clock / (double)Creatures[Percepteur].InternalClock);
            Creatures[Percepteur].Inputs[55] = TassementEntree(Creatures[Percepteur].Alea);
            Creatures[Percepteur].Inputs[56] = TassementEntree(Creatures[Percepteur].TempsDepuisDerniereReproduction);
            Creatures[Percepteur].Inputs[57] = TassementEntree(Creatures[Percepteur].TempsDepuisDernierClonage);
        }
        void CreateTerrainInputs(int Percepteur)
        {
            double x_ = Creatures[Percepteur].Position.Real;
            double y_ = Creatures[Percepteur].Position.Imaginary;
            Complex K = GetXYofTileAt(new Complex(x_, y_));
            //A l aide de la fonction GetFoccus
            //Terrain 
            Complex Foccus;
            Complex Viande = 0;
            Complex Herbe = 0;
            Complex Encombrement = 0;
            Complex Altitude = 0;
            Complex Coloration = 0;
            Complex Nutritif = 0;
            //On parcourt un cercle autour de la creature
            int CasesY = (int)Math.Round(CE.DistanceBrouillard / Scale);
            int Y0 = (int)Math.Max(K.Imaginary - CasesY,0);
            int Y1 = (int)Math.Min(K.Imaginary + CasesY,TerrainSubdiv-1);
            for (int y = Y0; y <= Y1; y++)
            {
                double u = Math.Pow(CasesY / (double)TerrainSubdiv, 2);
                double v = Math.Pow(Math.Abs(K.Imaginary - y) / (double)TerrainSubdiv, 2);
                double delta = u - v;
                int CasesX=(int)(TerrainSubdiv*Math.Sqrt(delta));
                int X0 = (int)Math.Max(K.Real - CasesX, 0);
                int X1 = (int)Math.Min(K.Real + CasesX, TerrainSubdiv - 1);
                for (int x = X0; x <= X1; x++)
                {
                    Complex PosTile = GetPositionOfTile(x, y);
                    Foccus = GetFoccus(Percepteur, PosTile) * (1.0 - Obstruction(Percepteur, PosTile)) * AccuiteVisuelle(Percepteur);
                    Viande += Foccus * Terrain[x, y].Viande;
                    Herbe += Foccus * Terrain[x, y].Vegetation;
                    Encombrement += Foccus * Terrain[x, y].CapaciteBloquage;
                    Altitude += Foccus * Terrain[x, y].Altitude;
                    Coloration += Foccus * Terrain[x, y].Coloration;
                    Nutritif += Foccus * Terrain[x, y].VegetationPouvoirNutritif;
                }
            }
            //On ramene en terme de densite sur la case
            Viande = Viande / (Scale * Scale);
            Herbe = Herbe / (Scale * Scale);
            Encombrement = Encombrement / (Scale * Scale);
            Creatures[Percepteur].Inputs[23] = TassementEntree(Viande);
            Creatures[Percepteur].Inputs[24] = TassementEntree(Herbe);
            Creatures[Percepteur].Inputs[25] = TassementEntree(Encombrement);
            Creatures[Percepteur].Inputs[52] = TassementEntree(Altitude);
            Creatures[Percepteur].Inputs[54] = TassementEntree(Coloration);
            Creatures[Percepteur].Inputs[62] = TassementEntree(Nutritif);
            //attention , si ajout ou retroiat d un input, penser a changer la valeur de NombreInputs dans constantes biologiques
        }
        void CreateRelativeInputs(int Percepteur)
        {
            //Creation des variables
            Complex Foccus;
            //Autres Individus
            //1 si peut se reproduire avec indiv, 0 sinon;
            Complex Reproductibilite = 0;
            //Autres variables relatives
            Complex SanteRelative = 0;
            Complex EnergieRelative = 0;Complex VitesseRelative = 0;
            Complex VitesseRotationRelative=0;
            Complex OrientationRelative = 0;
            //Renvoie 1 si individu present
            Complex Densite = 0;Complex AttitudeRelative = 0;
            Complex VolonteReproductiveRelative = 0;Complex AgresiviteRelative = 0;
            Complex FoccusLineaireRelatif = 0;Complex FoccusAngulaireRelatif = 0;
            //Variables determinees
            Complex TempsDepuisDerniereReproductionRelatif = 0;
            Complex TempsDepuisDernierClonageRelatif = 0;
            Complex AgeRelatif = 0;Complex MasseRelative = 0;
            Complex TemperaturePercueRelative = 0;Complex EffortRelatif = 0;
            //Stats
            Complex DefenseRelative = 0;Complex AttaqueRelative = 0;
            Complex MaxVitesseRelative = 0;Complex CouleurRelative = 0;
            Complex MaxVitesseRotationRelative = 0;
            Complex RegimeRelatif = 0;Complex GuerisonRelative = 0;
            Complex IdealTempRelative = 0;Complex IdealLumiereRelative = 0;
            Complex DonNRJEnfantRelatif = 0;
            //Calcul des variables
            for (int i = 0; i < Creatures.GetLongLength(0); i++)
            {
                
                    if (EnVie[i] && i != Percepteur)
                    {
                        Complex pos = Creatures[i].Position;
                        Foccus = GetFoccus(Percepteur, pos) * (1.0-Obstruction(Percepteur, pos)) * (1.0-Camouflage(i)) * AccuiteVisuelle(Percepteur);
                        //Variables particulieres
                        Densite += Foccus;
                        if (Creatures[Percepteur].Brain.Compatible(Creatures[i].Brain))
                        {
                            Reproductibilite += Foccus;
                        }
                        VitesseRelative += Foccus * Math.Abs((Creatures[Percepteur].Vitesse.Magnitude - Creatures[i].Vitesse.Magnitude));
                        OrientationRelative += Foccus * Math.Abs((Creatures[Percepteur].Vitesse.Phase - Creatures[i].Vitesse.Phase));

                    //Variables classiques
                        VitesseRotationRelative += Foccus * (Creatures[Percepteur].VitesseRotation - Creatures[i].VitesseRotation);
                        SanteRelative += Foccus * (Creatures[Percepteur].Sante - Creatures[i].Sante);
                        EnergieRelative += Foccus * (Creatures[Percepteur].Energie - Creatures[i].Energie);
                        AttitudeRelative += Foccus * (Creatures[Percepteur].Attitude - Creatures[i].Attitude);
                        VolonteReproductiveRelative += Foccus * (Creatures[Percepteur].VolonteReproductive - Creatures[i].VolonteReproductive);
                        AgresiviteRelative += Foccus * (Creatures[Percepteur].Agresivite - Creatures[i].Agresivite);
                        FoccusLineaireRelatif += Foccus * (Creatures[Percepteur].FoccusLineaire - Creatures[i].FoccusLineaire);
                        FoccusAngulaireRelatif += Foccus * (Creatures[Percepteur].FoccusAngulaire - Creatures[i].FoccusAngulaire);
                        TempsDepuisDerniereReproductionRelatif += Foccus * (Creatures[Percepteur].TempsDepuisDerniereReproduction - Creatures[i].TempsDepuisDerniereReproduction);
                        TempsDepuisDernierClonageRelatif += Foccus * (Creatures[Percepteur].TempsDepuisDernierClonage - Creatures[i].TempsDepuisDernierClonage);
                        AgeRelatif += Foccus * (Creatures[Percepteur].Age - Creatures[i].Age);
                        MasseRelative += Foccus * (Creatures[Percepteur].Masse - Creatures[i].Masse);
                        TemperaturePercueRelative += Foccus * (Creatures[Percepteur].TemperaturePercue - Creatures[i].TemperaturePercue);
                        EffortRelatif += Foccus * (Creatures[Percepteur].Effort - Creatures[i].Effort);
                        DefenseRelative += Foccus * (Creatures[Percepteur].Defense - Creatures[i].Defense);
                        AttaqueRelative += Foccus * (Creatures[Percepteur].Attaque - Creatures[i].Attaque);
                        MaxVitesseRelative += Foccus * (Creatures[Percepteur].MaxVitesse - Creatures[i].MaxVitesse);
                        MaxVitesseRotationRelative += Foccus * (Creatures[Percepteur].MaxVitesseRotation - Creatures[i].MaxVitesseRotation);
                        CouleurRelative += Foccus * (Creatures[Percepteur].Couleur - Creatures[i].Couleur);
                        RegimeRelatif += Foccus * (Creatures[Percepteur].Regime - Creatures[i].Regime);
                        GuerisonRelative += Foccus * (Creatures[Percepteur].Guerison - Creatures[i].Guerison);
                        IdealTempRelative += Foccus * (Creatures[Percepteur].IdealTemp - Creatures[i].IdealTemp);
                        IdealLumiereRelative += Foccus * (Creatures[Percepteur].LumiereIdeale - Creatures[i].LumiereIdeale);
                        DonNRJEnfantRelatif += Foccus * (Creatures[Percepteur].DonNRJEnfant - Creatures[i].DonNRJEnfant);
                    }
                
            }
            //Attribution des valeurs
            Creatures[Percepteur].Inputs[26] = TassementEntree(SanteRelative);
            Creatures[Percepteur].Inputs[27] = TassementEntree(EnergieRelative);
            Creatures[Percepteur].Inputs[28] = TassementEntree(VitesseRelative);
            Creatures[Percepteur].Inputs[29] = TassementEntree(Densite);
            Creatures[Percepteur].Inputs[30] = TassementEntree(AttitudeRelative);
            Creatures[Percepteur].Inputs[31] = TassementEntree(VolonteReproductiveRelative);
            Creatures[Percepteur].Inputs[32] = TassementEntree(AgresiviteRelative);
            Creatures[Percepteur].Inputs[33] = TassementEntree(FoccusAngulaireRelatif);
            Creatures[Percepteur].Inputs[34] = TassementEntree(FoccusLineaireRelatif);
            Creatures[Percepteur].Inputs[35] = TassementEntree(TempsDepuisDerniereReproductionRelatif);
            Creatures[Percepteur].Inputs[36] = TassementEntree(AgeRelatif);
            Creatures[Percepteur].Inputs[37] = TassementEntree(MasseRelative);
            Creatures[Percepteur].Inputs[38] = TassementEntree(TemperaturePercueRelative);
            Creatures[Percepteur].Inputs[39] = TassementEntree(EffortRelatif);
            Creatures[Percepteur].Inputs[40] = TassementEntree(DefenseRelative);
            Creatures[Percepteur].Inputs[41] = TassementEntree(AttaqueRelative);
            Creatures[Percepteur].Inputs[42] = TassementEntree(MaxVitesseRelative);
            Creatures[Percepteur].Inputs[43] = TassementEntree(CouleurRelative);
            Creatures[Percepteur].Inputs[44] = TassementEntree(RegimeRelatif);
            Creatures[Percepteur].Inputs[45] = TassementEntree(GuerisonRelative);
            Creatures[Percepteur].Inputs[46] = TassementEntree(IdealTempRelative);
            Creatures[Percepteur].Inputs[47] = TassementEntree(DonNRJEnfantRelatif);
            Creatures[Percepteur].Inputs[49] = TassementEntree(Reproductibilite);
            Creatures[Percepteur].Inputs[50] = TassementEntree(IdealLumiereRelative);
            Creatures[Percepteur].Inputs[58] = TassementEntree(TempsDepuisDernierClonageRelatif);
            Creatures[Percepteur].Inputs[59] = TassementEntree(MaxVitesseRotationRelative);
            Creatures[Percepteur].Inputs[60] = TassementEntree(OrientationRelative);
            Creatures[Percepteur].Inputs[61] = TassementEntree(VitesseRotationRelative);
        }
        void RunNetwork(int Percepteur)
        {
            Creatures[Percepteur].Brain.EffectuerCalcul(Creatures[Percepteur].Inputs, ref Creatures[Percepteur].Outputs);

        }
        void CreateInputs(int Percepteur)
        {
            CreateInternalInputs(Percepteur);
            CreateTerrainInputs(Percepteur);
            CreateRelativeInputs(Percepteur);
            for(int i=0; i<CB.NombreInputs;i++)
            {
                double activation;
                if(Creatures[Percepteur].Perception[i]>=1.0)
                {
                    activation = 1.0;
                }
                else
                {
                    activation = 0.0;
                }
                Creatures[Percepteur].Inputs[i] = Creatures[Percepteur].Inputs[i] * activation;
            }
            // 63 inputs
        }
        //Update Princiapux
        void UpdatePerception(int Percepteur)
        {
            UpdateVision(Percepteur);
            UpdateVentPercu(Percepteur);
            UpdateTemperaturePercue(Percepteur);
            UpdateHorloges(Percepteur);
        }
        void UpdateVariablesInternes(int Percepteur)
        {
            UpdateMasse(Percepteur);
            UpdateRayon(Percepteur);
            UpdateEffort(Percepteur);
        }
        void UpdateAction(int Percepteur, ref Random r)
        {
            GenererAlea(Percepteur, ref r);
            UpdateAlimentation(Percepteur);
            UpdateAgressivite(Percepteur);
            UpdateAttitude(Percepteur);
            UpdateVolonteReproductive(Percepteur);
            UpdateVitesse(Percepteur);
            Mouvement(Percepteur);
            SanteUpdateEtCout(Percepteur);
            AppliquerBilanEnergetique(Percepteur);
            CombatsEtReproduction(Percepteur, ref r);
        }
        void UpdateDepenses(int Percepteur)
        {
            TemperatureDepenses(Percepteur);
            MouvementDepense(Percepteur);
            SubsistanceDepense(Percepteur);
        }
        void UpdateTerrain(ref Random r)
        {
            VE.Time+= dt;
            //La lumiere a pour periode la duree d'un jour et est comprise entre 0 et 2;
            VE.Lumiere = 1.0+Math.Cos((float)VE.Time*(2.0*Math.PI)/CE.DureeJour);
            //Le vent est modifié chaque frame
            Complex dVent = Complex.FromPolarCoordinates(CE.VariabiliteVent* r.NextDouble(), r.NextDouble()*2.0*Math.PI);
            VE.Vent += dVent*dt;
            //La temperature varie proportionellement au vent et aussi d'une composante aleatoire
            VE.Temperature += dt *((VE.Lumiere-1.0)*VE.Vent.Magnitude+CE.VariabiliteTemp*(r.NextDouble()-0.5));
            VE.Temperature = Math.Max(0.5, VE.Temperature);
            for (int x = 0; x < Terrain.GetLength(0); x++)
            {
                for (int y = 0; y < Terrain.GetLength(1); y++)
                {
                    //On fait pousser les plantes et pourrir la viande
                    double VGIdealT = Terrain[x, y].IdealTempVegetation;
                    double TempVGFactor = (Math.Exp(-1.0 * (Math.Pow(VE.Temperature - VGIdealT, 2) * Terrain[x, y].VegetationTempSensibilite)));
                    double Cvg = dt * CE.CroissanceVegetaux * (VE.Lumiere * Terrain[x, y].Fertilité * TempVGFactor-CE.SeuilMortVegetaux);
                    Terrain[x, y].Vegetation += Cvg; 
                    Terrain[x, y].Vegetation =Math.Max(Terrain[x, y].Vegetation,0);
                    double Viande = Terrain[x, y].Viande;
                    double VPerte = Viande * CE.TauxDegradationViande*(1.0 - Math.Exp(-CE.Pourriture * VE.Temperature))*dt;
                    Terrain[x, y].Viande = Math.Max(Terrain[x, y].Viande -VPerte,0);
                    //L'encombrement est proportionel a la vegetation et l 'impraticabilite
                    double v = Terrain[x, y].Vegetation * Terrain[x, y].VegetationEncombrement;
                    double I = Terrain[x, y].Impraticabilité;
                    Terrain[x, y].CapaciteBloquage =  v+I ;
                    //Couleur apparente du terrain
                    double CoeffHerbe = (1.0-1.0/(1.0+ Terrain[x, y].Vegetation*CE.ImportanceColorationVegetaux/(Scale*Scale)));
                    Terrain[x, y].Coloration = (Terrain[x, y].CouleurTerrainNu+CoeffHerbe* Terrain[x, y].CouleurVegetation)/(1.0+CoeffHerbe);
                }
            }
        }
        void CreateInputsAndRunNetwork(Object InfoP)
        {
            IndexInfo Info = (IndexInfo)InfoP;
            for (int i = Info.Start; i < Info.Stop; i++)
            {
                if (EnVie[i])
                {
                    CreateInputs(i);
                    RunNetwork(i);
                }
            }
            
        }
        //Update general
        public void TourSuivant (Creature MaxValues)
        {
            log = "";

            PerteThermiques = 0;
            PertesMouvement = 0;
            PertesSubsistance = 0;
            PertesSante = 0;
            PertesReproduction = 0;
            PertesCombat = 0;
            GainAlimentation = 0;
            System.Diagnostics.Stopwatch s = new System.Diagnostics.Stopwatch();
            //Updates Generales
            //On ajoute des individus aleatoires
            AjoutIndividuAleatoire(ProbaPoisson(GP.FrequenceApparitionNouvelIndividu),ref GP.RANDOM,MaxValues);
            s.Start();
            //Le terrain est mis a jour
            UpdateTerrain(ref GP.RANDOM);
            s.Stop();
            log+=("---------" + "Update Terrain " + s.Elapsed+"\n");
            s.Reset();
            s.Start();
            //Variables Internes
            for (int i=0; i<Creatures.GetLength(0); i++)
            {
                if (EnVie[i])
                {
                    UpdateVariablesInternes(i);
                    UpdatePerception(i);
                    UpdateDepenses(i);
                }
            }
            IndexInfo I1 = new IndexInfo();
            IndexInfo I2 = new IndexInfo();
            IndexInfo I3 = new IndexInfo();
            IndexInfo I4 = new IndexInfo();
            s.Stop();
            log+=("---------" + "Variables Internes et depenses " + s.Elapsed+"\n");
            s.Reset();
            s.Start();
           for (int i = 0; i < Creatures.GetLength(0); i++)
            {
                if (EnVie[i])
                {
                    UpdateAction(i, ref GP.RANDOM);
                }
            }
            s.Stop();
            log+=("---------" + "Actions " + s.Elapsed+"\n");
            s.Reset();
            s.Start();
            int Middle = Creatures.GetLength(0)/2;
            int quarter = Middle / 2;
            if (Creatures.GetLength(0) > 10)
            {
                //MUltithread
                I1.Start = 0;
                I1.Stop = quarter;
                I2.Start = quarter;
                I2.Stop = Middle;
                I3.Start = Middle;
                I3.Stop = Middle + quarter;
                I4.Start = Middle + quarter;
                I4.Stop = Creatures.GetLength(0);


                //On lance en parallele les inputs et le network
                Thread Cmod0 = new Thread(new ParameterizedThreadStart(CreateInputsAndRunNetwork));
                Thread Cmod1 = new Thread(new ParameterizedThreadStart(CreateInputsAndRunNetwork));
                Thread Cmod2 = new Thread(new ParameterizedThreadStart(CreateInputsAndRunNetwork));
                Thread Cmod3 = new Thread(new ParameterizedThreadStart(CreateInputsAndRunNetwork));
                Cmod0.Start(I1);
                Cmod1.Start(I2);
                Cmod2.Start(I3);
                Cmod3.Start(I4);
                Cmod3.Join();
                Cmod2.Join();
                Cmod1.Join();
                Cmod0.Join();
            }
            else
            {
                //Monothread
                I1.Start = 0;
                I1.Stop = Creatures.GetLength(0);
                CreateInputsAndRunNetwork(I1);
            }
            //ICI
            s.Stop();
            log+=("---------" + "Variables Neuronales " + s.Elapsed+"\n");
            s.Reset();
            s.Start();
            PopulationUpdate();
            s.Stop();
            log+=("---------" + "Population Update " + s.Elapsed+"\n");
            Enfants = new List<Creature>();
            double N = Creatures.Count();
            double ptmoy = PerteThermiques / N;
            double pmovmoy = PertesMouvement / N;
            double psubmoy = PertesSubsistance / N;
            double psamoy = PertesSante / N;
            double prepmoy = PertesReproduction / N;
            double pcombmoy = PertesCombat / N;
            double galmoy = GainAlimentation / N;
            double bilan = galmoy - pcombmoy - prepmoy - psamoy - psubmoy - pmovmoy - ptmoy;
            log+=(">>>>>>>Energies Moyennes :"+"\n");
            log+=("Thermique : "+ptmoy + "\n");
            log+=("Mouvement : "+pmovmoy + "\n");
            log+=("Subsistance : "+psubmoy + "\n");
            log+=("Santé : " + psamoy + "\n");
            log+=("Reproduction : " + prepmoy + "\n");
            log+=("Combat : " + pcombmoy + "\n");
            log+=("Alimentation : " + galmoy + "\n");
            log+=("Bilan : " + bilan + "\n");
            
            
        }
        //Generation
        void GenerateTerrain(int X, int Y, ref Random R, ElementTerrain MaxTerrainValues)
        {
            Terrain = new ElementTerrain[X,Y]; 
            float Frequency = 0.02f;
            float Amplitude = 0.8f;
            float Persistance = 0.7f;
            int Octaves = 7;
            Perlin perlin1 = new Perlin(X, Y, ref R);
            Perlin perlin2 = new Perlin(X, Y, ref R);
            Perlin perlin3 = new Perlin(X, Y, ref R);
            Perlin perlin4 = new Perlin(X, Y, ref R);
            Perlin perlin5 = new Perlin(X, Y, ref R);
            Perlin perlin6 = new Perlin(X, Y, ref R);
            Perlin perlin7 = new Perlin(X, Y, ref R);
            Perlin perlin8 = new Perlin(X, Y, ref R);
            Perlin perlin9 = new Perlin(X, Y, ref R);
            Perlin perlin10 = new Perlin(X, Y, ref R);
            for (int x = 0; x < TerrainSubdiv; x++)
            {
                for (int y = 0; y < TerrainSubdiv; y++)
                {
                    double P1 = perlin1.GenerateValues(x, y, (float)MaxTerrainValues.Fertilité, Frequency, Amplitude, Persistance, Octaves);
                    double P2 = perlin2.GenerateValues(x, y, (float)MaxTerrainValues.IdealTempVegetation, Frequency, Amplitude, Persistance, Octaves);
                    double P3 = perlin3.GenerateValues(x, y, (float)MaxTerrainValues.VegetationTempSensibilite, Frequency, Amplitude, Persistance, Octaves);
                    double P4 = perlin4.GenerateValues(x, y, (float)MaxTerrainValues.Impraticabilité, Frequency, Amplitude, Persistance, Octaves);
                    double P5 = perlin5.GenerateValues(x, y, (float)MaxTerrainValues.Vegetation, Frequency, Amplitude, Persistance, Octaves);
                    double P6 = perlin6.GenerateValues(x, y, (float)MaxTerrainValues.Altitude, Frequency, Amplitude, Persistance, Octaves);
                    double P7 = perlin7.GenerateValues(x, y, (float)MaxTerrainValues.CouleurTerrainNu, Frequency, Amplitude, Persistance, Octaves);
                    double P8 = perlin8.GenerateValues(x, y, (float)MaxTerrainValues.CouleurVegetation, Frequency, Amplitude, Persistance, Octaves);
                    double P9 = perlin9.GenerateValues(x, y, (float)MaxTerrainValues.VegetationEncombrement, Frequency, Amplitude, Persistance, Octaves);
                    double P10 = perlin10.GenerateValues(x, y, (float)MaxTerrainValues.VegetationPouvoirNutritif, Frequency, Amplitude, Persistance, Octaves);
                    Terrain[x, y].Fertilité = P1;
                    Terrain[x, y].IdealTempVegetation = P2;
                    Terrain[x, y].VegetationTempSensibilite = P3;
                    Terrain[x, y].Impraticabilité = P4;
                    Terrain[x, y].Vegetation = P5;
                    Terrain[x, y].Altitude = P6;
                    Terrain[x, y].CouleurTerrainNu = P7;
                    Terrain[x, y].CouleurVegetation = P8;
                    Terrain[x, y].VegetationEncombrement = P9;
                    Terrain[x, y].VegetationPouvoirNutritif = P10;
                }
            }
        }
        void GenererPremiersPeuples(int NombreIndividus, ref Random r, Creature MaxValues)
        {
            Creatures = new Creature[NombreIndividus];
            EnVie = new bool[NombreIndividus];
            VE.IDlast = 0;
            for (int i = 0; i < NombreIndividus; i++)
            {
                
                EnVie[i] = true;
                Creature c = new Creature();
                GenererCreatureAleatoire(ref c,ref r,MaxValues);
                Creatures[i] = c;
            }
        }
        //ValeursMax definies ici
        void GenererCreatureAleatoire(ref Creature C,ref Random r, Creature MaxValues)
        {

            VE.IDlast++;

            C = new Creature();
            C.Sante = r.NextDouble();
            C.Position = new Complex(r.NextDouble(), r.NextDouble());
            C.Energie = r.NextDouble()*MaxValues.Energie;
            C.TempsDepuisDerniereReproduction = 0;
            C.TempsDepuisDernierClonage = 0;
            C.Alea = 0;
            C.Age = 0;
            C.Couleur = r.NextDouble();
            C.P1Id = 0;
            C.P2Id = 0;
            C.nom = new Nom();
            C.nom.Generer(ref r);
            C.ID = VE.IDlast;
            C.Regime = Math.Max(0,Math.Min(1,0.5 + (r.NextDouble() - 0.5) * MaxValues.Regime));
            C.MaxVitesse = MaxValues.MaxVitesse * r.NextDouble();
            C.Vitesse = Complex.FromPolarCoordinates(r.NextDouble()*C.MaxVitesse, Math.PI * 2 * r.NextDouble());
            C.MaxVitesseRotation = MaxValues.MaxVitesseRotation * r.NextDouble();
            C.InternalClock = Math.Max(Convert.ToInt32(r.NextDouble() * MaxValues.InternalClock), 1);
            C.Guerison = MaxValues.Guerison * (r.NextDouble() - 0.5);
            C.DonNRJEnfant = r.NextDouble() * MaxValues.DonNRJEnfant;
            C.Attaque = r.NextDouble() * MaxValues.Attaque;
            C.Defense = r.NextDouble() * MaxValues.Defense;
            C.ReproducedThisTurn = false;
            C.IdealTemp = r.NextDouble() * MaxValues.IdealTemp;
            C.LumiereIdeale = Math.Min(0,Math.Max(0, 1.0 + (2.0 * r.NextDouble() - 1) * MaxValues.LumiereIdeale));
            int nombreInputs = CB.NombreInputs;
            int nombreOutputs = CB.NombreOutputs;
            C.Inputs = new Complex[nombreInputs];
            for(int z=0; z< nombreInputs; z++)
            {
                C.Inputs[z] = Complex.FromPolarCoordinates(r.NextDouble(),Math.PI*2*r.NextDouble());
            }
            C.Outputs = new Complex [nombreOutputs];
            for (int z = 0; z < nombreOutputs; z++)
            {
                C.Outputs[z] = Complex.FromPolarCoordinates(r.NextDouble(), Math.PI * 2 * r.NextDouble());
            }
            //Construire Cerveau
            int Layers = CouchesDuReseauEtVariabilite.GetLength(0);
            int[] Struct = new int[Layers+2];
            //Inputs
            Struct[0] = nombreInputs;
            
            //Couches Intermediaires
            for (int i = 0; i < Layers; i++)
            {
                Struct[i + 1] = CouchesDuReseauEtVariabilite[i, 0] + r.Next(0, CouchesDuReseauEtVariabilite[i, 1]);
            }
            //Outputs
            Struct[Layers + 1] = nombreOutputs;
            Cerveau Net = new Cerveau();
            Net.GenererAleatoire(Struct, ref r);
            C.Brain = Net;
            //Organes de perception
            C.Perception = new double[nombreInputs];
            for (int i = 0; i < nombreInputs; i++)
            {
                C.Perception[i] = 1+GaussianRandom(ref r)*GP.DiversiteLuciditeInitiale;
            }

        }
        Creature ReproductionSexuee(int IndivA, int IndivB,double TauxMutation,  ref Random r)
        {
            //L enfant est inspire de ses parents et d'une variable aleatoire
            Creature Enfant = new Creature();
            VE.IDlast++;
            Enfant.Age = 0;
            Enfant.Sante = 1.0;
            Enfant.P1Id = Creatures[IndivA].ID;
            Enfant.P2Id = Creatures[IndivB].ID;
            Enfant.ID = VE.IDlast;
            
            double j = (Creatures[IndivA].DonNRJEnfant) * Creatures[IndivA].Energie;
            double k = (Creatures[IndivB].DonNRJEnfant) * Creatures[IndivB].Energie;
            Enfant.Energie = (j+k) * CB.RendementReproductionSexuee;
            Creatures[IndivA].BilanEnergetique -= j;
            Creatures[IndivB].BilanEnergetique -= k;
            PertesReproduction += j + k;
            Enfant.Position = 0.5 * (Creatures[IndivA].Position + Creatures[IndivB].Position);
            Enfant.Vitesse = 0;
            Enfant.nom = new Nom();
            Enfant.nom.Croiser(ref r, Creatures[IndivA].nom, Creatures[IndivB].nom,GP.ProbaChangementNomDeFamille);
            log+=(">"+Enfant.nom.Prenom + " " + Enfant.nom.NomDeFamille + " est né"+"\n");
            Enfant.Defense = Math.Max(0, GaussianRandom(ref r)*TauxMutation+ (Creatures[IndivA].Defense + Creatures[IndivB].Defense)/2);
            Enfant.Attaque = Math.Max(0, GaussianRandom(ref r)*TauxMutation + (Creatures[IndivA].Attaque + Creatures[IndivB].Attaque) / 2);
            Enfant.MaxVitesse = Math.Max(0, GaussianRandom(ref r)*TauxMutation + (Creatures[IndivA].MaxVitesse + Creatures[IndivB].MaxVitesse) / 2);
            Enfant.MaxVitesseRotation = Math.Max(0, GaussianRandom(ref r) * TauxMutation + (Creatures[IndivA].MaxVitesseRotation + Creatures[IndivB].MaxVitesseRotation) / 2);
            Enfant.Couleur = (Creatures[IndivA].Couleur + Creatures[IndivB].Couleur+r.NextDouble()*TauxMutation)/(2+TauxMutation) ;
            Enfant.Regime = (Creatures[IndivA].Regime + Creatures[IndivB].Regime + r.NextDouble() * TauxMutation) / (2 + TauxMutation);
            Enfant.DonNRJEnfant = (Creatures[IndivA].DonNRJEnfant + Creatures[IndivB].DonNRJEnfant + r.NextDouble()*TauxMutation) / (2+TauxMutation);
            Enfant.Guerison =  GaussianRandom(ref r)*TauxMutation + (Creatures[IndivA].Guerison + Creatures[IndivB].Guerison)/2.0;
            Enfant.IdealTemp = Math.Max(0, GaussianRandom(ref r)*TauxMutation + (Creatures[IndivA].IdealTemp + Creatures[IndivB].IdealTemp) / 2);
            Enfant.LumiereIdeale = Math.Min(2,Math.Max(0,(Creatures[IndivA].LumiereIdeale + Creatures[IndivB].LumiereIdeale + 2.0*r.NextDouble()*TauxMutation) / (2+TauxMutation)));
            //Prend le modele des parents pour les in outputs
            Enfant.Inputs = Creatures[IndivA].Inputs;
            Enfant.Outputs = Creatures[IndivA].Outputs;
            Enfant.InternalClock = (int)Math.Max(1, GaussianRandom(ref r)*TauxMutation + (Creatures[IndivA].InternalClock + Creatures[IndivB].InternalClock) / 2);
            Enfant.Brain = new Cerveau();
            Enfant.Brain = Enfant.Brain.Moyenne(Creatures[IndivA].Brain,Creatures[IndivB].Brain);
            Enfant.Brain.Mutations(GP.TauxMutationReproductionSexuee, ref r);
            Enfant.Perception = new double[CB.NombreInputs];
            for (int i = 0; i < CB.NombreInputs; i++)
            {
                Enfant.Perception[i] = GaussianRandom(ref r) * TauxMutation / CB.StabiliteOrganesPerception+(Creatures[IndivA].Perception[i]+ Creatures[IndivB].Perception[i])/2.0;
            }
            return Enfant;

        }
        Creature ReproductionClonage(int Indiv, double TauxMutation, ref Random r)
        {
            Creature Enfant = new Creature();
            VE.IDlast++;
            Enfant.Age = 0;
            Enfant.Sante = 1.0;
            Enfant.P1Id = Creatures[Indiv].ID;
            Enfant.P2Id = Creatures[Indiv].ID;
            Enfant.ID = VE.IDlast;
            AppliquerBilanEnergetique(Indiv);
            double EnergieDisp = Creatures[Indiv].Energie * CB.RendementClonage;
            Enfant.Energie = 0.5 *EnergieDisp ;
            PertesReproduction += Creatures[Indiv].Energie * (1.0 - 0.5 * CB.RendementClonage);
            Creatures[Indiv].Energie = 0.5 * EnergieDisp;
            Enfant.Position = Creatures[Indiv].Position;
            Enfant.Vitesse = 0;
            Enfant.nom = new Nom();
            Enfant.nom.Croiser(ref r, Creatures[Indiv].nom, Creatures[Indiv].nom, GP.ProbaChangementNomDeFamille);
            log+=(">" + Enfant.nom.Prenom + " " + Enfant.nom.NomDeFamille + "est né par clonage de " + Creatures[Indiv].nom.Prenom + " " + Creatures[Indiv].nom.NomDeFamille+"\n");
            Enfant.Defense = Math.Max(0, GaussianRandom(ref r) * TauxMutation + Creatures[Indiv].Defense);
            Enfant.Attaque = Math.Max(0, GaussianRandom(ref r) * TauxMutation + Creatures[Indiv].Attaque);
            Enfant.MaxVitesse = Math.Max(0, GaussianRandom(ref r) * TauxMutation + Creatures[Indiv].MaxVitesse);
            Enfant.MaxVitesseRotation = Math.Max(0, GaussianRandom(ref r) * TauxMutation + Creatures[Indiv].MaxVitesseRotation);
            Enfant.Couleur = (Creatures[Indiv].Couleur + r.NextDouble() * TauxMutation) / (1 + TauxMutation);
            Enfant.Regime = (Creatures[Indiv].Regime + r.NextDouble() * TauxMutation) / (1 + TauxMutation);
            Enfant.DonNRJEnfant = (Creatures[Indiv].DonNRJEnfant + r.NextDouble() * TauxMutation) / (1 + TauxMutation);
            Enfant.Guerison = GaussianRandom(ref r) * TauxMutation + Creatures[Indiv].Guerison;
            Enfant.IdealTemp = Math.Max(0, GaussianRandom(ref r) * TauxMutation + Creatures[Indiv].IdealTemp);
            Enfant.LumiereIdeale = Math.Min(2, Math.Max(0, (Creatures[Indiv].LumiereIdeale + 2.0 * r.NextDouble() * TauxMutation) / (1 + TauxMutation)));
            //Prend le modele des parents pour les in outputs
            Enfant.Inputs = Creatures[Indiv].Inputs;
            Enfant.Outputs = Creatures[Indiv].Outputs;
            Enfant.InternalClock = (int)Math.Max(1, GaussianRandom(ref r) * TauxMutation + Creatures[Indiv].InternalClock);
            Enfant.Brain = new Cerveau();
            Enfant.Brain = Creatures[Indiv].Brain;
            Enfant.Brain.Mutations(GP.TauxMutationClonage, ref r);
            Enfant.Perception = Creatures[Indiv].Perception;
            for(int i=0;i<CB.NombreInputs;i++)
            {
                Enfant.Perception[i] += GaussianRandom(ref r) * TauxMutation/CB.StabiliteOrganesPerception;
            }
            return Enfant;
        }
        //Evenements
        void GenererAlea(int Percepteur, ref Random r)
        {
            Creatures[Percepteur].Alea = Complex.FromPolarCoordinates(r.NextDouble(), 2 * Math.PI * r.NextDouble());
        }
        //Les evenements ont ils lieu, et resolution
        bool Combat(int IndivA, int IndivB, ref Random r)
         {
             bool CombatALieu = false;
             double BlessureB = 0;
            //Si l'individu veut Attaquer
            double AgiliteAttaque = Creatures[IndivA].Vitesse.Magnitude + Creatures[IndivA].VitesseRotation;
            double AgiliteDefense = Creatures[IndivB].Vitesse.Magnitude + Creatures[IndivB].VitesseRotation;
            double Freq = Creatures[IndivA].Agresivite*GP.FrequenceMaxCombat*Sig(AgiliteAttaque-AgiliteDefense);
                if (r.NextDouble()<ProbaPoisson(Freq))
                {
                    CombatALieu = true;
                    //L'energie impliquée dans l'attaque
                    double EnergieEmployeeAttaque = Creatures[IndivA].Agresivite * Creatures[IndivA].Attaque * CB.CoutAttaque* Creatures[IndivA].Masse;
                    double EnergieEmployeeDefense = (1.0-Creatures[IndivB].Agresivite) * Creatures[IndivB].Defense * CB.CoutAttaque * Creatures[IndivB].Masse;
                //La blessure infligée dépend de la sante des individus
                BlessureB = Sig((EnergieEmployeeAttaque * Creatures[IndivA].Sante - EnergieEmployeeDefense * Creatures[IndivB].Sante));
                    //L'energie est depensee
                    Creatures[IndivA].BilanEnergetique -= EnergieEmployeeAttaque;
                    Creatures[IndivB].BilanEnergetique -= EnergieEmployeeDefense;
                PertesCombat += EnergieEmployeeAttaque + EnergieEmployeeDefense;
                    //On applique la blessure
                    Creatures[IndivB].Sante -= BlessureB;
                }
                return CombatALieu;
         }
        bool Reproduction(int IndivA, int IndivB, ref Random r)
        {
            //les creatures peuvent elles se reproduire (pas de reproduction a eu lieu ce tour ci) et sont elles compatibles neuronalement
            if ((Creatures[IndivA].Brain.Compatible(Creatures[IndivB].Brain))&&(!Creatures[IndivA].ReproducedThisTurn && !Creatures[IndivB].ReproducedThisTurn))
            {
                //L'importance de la decision est attribuee au prorata de la force
                double DominationA = Creatures[IndivA].Sante * (Creatures[IndivA].Attaque + Creatures[IndivA].Defense)* Creatures[IndivA].Masse;
                double DominationB = Creatures[IndivB].Sante * (Creatures[IndivB].Attaque + Creatures[IndivB].Defense) * Creatures[IndivB].Masse;
                //La probabilite depend de la volonte des individus (0.5 = indifferent)
                double Freq = GP.FrequenceMaxReproduction*Sig((Creatures[IndivA].VolonteReproductive * DominationA+Creatures[IndivB].VolonteReproductive*DominationB)/(DominationA+DominationB) );
                if (r.NextDouble() < ProbaPoisson(Freq))
                {
                    //La reproduction a lieu
                    Creatures[IndivA].ReproducedThisTurn = true;
                    Creatures[IndivB].ReproducedThisTurn = true;
                    //On remets les horloges a zero
                    Creatures[IndivA].TempsDepuisDerniereReproduction = 0;
                    Creatures[IndivB].TempsDepuisDerniereReproduction = 0;
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
        bool Clonage(int Indiv, ref Random r)
        {
            double Freq = GP.FrequenceMaxClonage*Creatures[Indiv].Outputs[7].Magnitude;
            if(r.NextDouble()<ProbaPoisson(Freq))
            {
                Creatures[Indiv].TempsDepuisDernierClonage = 0;
                return true;
            }
            else
            {
                return false;
            }
        }
            //Application des evenements
        void CombatsEtReproduction(int Percepteur, ref Random r)
        {
            if (EnVie[Percepteur])
            {
                //Clonage
                if(Clonage(Percepteur,ref r))
                {
                    Enfants.Add(ReproductionClonage(Percepteur, GP.TauxMutationClonage, ref r));
                }
                //On cherche la plus proche des creatures en contact
                //Quelle est la creature la plus proche
                double MinDist = double.PositiveInfinity;
                double dist = 0;
                int MinIndex = 0;
                for (int i = 0; i < Creatures.GetLongLength(0); i++)
                {
                    if (EnVie[i] && i != Percepteur)
                    {
                        //Distance exprimee en metres
                        dist = Distance(Creatures[Percepteur].Position, Creatures[i].Position);//en m
                        if (dist < MinDist)
                        {
                            MinDist = dist;
                            MinIndex = i;
                        }
                    }
                }
                Creatures[Percepteur].SEstBattu = false;
                //Cette creature est elle ne contact ?
                if (MinDist < Math.Max(CE.RayonInteraction, Creatures[Percepteur].Rayon + Creatures[MinIndex].Rayon))
                {
                    Nom nomA = Creatures[Percepteur].nom;
                    Nom nomB = Creatures[MinIndex].nom;
                    if (EnVie[Percepteur] && EnVie[MinIndex])
                    {
                        //Contact
                        if (Reproduction(Percepteur, MinIndex, ref r))
                        {
                            //La reproduction a eu lieu
                            log+=("<>"+nomA.Prenom + " " + nomA.NomDeFamille + " se reproduit avec " + nomB.Prenom + " " + nomB.NomDeFamille+"\n");
                            Enfants.Add(ReproductionSexuee(Percepteur, MinIndex,GP.TauxMutationReproductionSexuee, ref r));
                        }
                    }
                    UpdateVivant(Percepteur, false);
                    UpdateVivant(MinIndex, false);
                    //On teste le combat potentiel
                    if (EnVie[Percepteur] && EnVie[MinIndex])
                    {
                        Creatures[Percepteur].SEstBattu = Combat(Percepteur, MinIndex, ref r);
                        if (Creatures[Percepteur].SEstBattu)
                        {
                            log+=("/"+nomA.Prenom + " " + nomA.NomDeFamille + " attaque " + nomB.Prenom + " " + nomB.NomDeFamille+"\n");
                        }
                    }
                    UpdateVivant(MinIndex, true);
                }
                if (EnVie[Percepteur])
                {
                    UpdateVivant(Percepteur, true);
                }
            }
        }
        //Mouvement Physique
        void Mouvement(int Percepteur)
        {
            //dx et dy sont exprimes en metres
            double dx = (Creatures[Percepteur].Vitesse.Real * dt );
            double dy = (Creatures[Percepteur].Vitesse.Imaginary * dt );
            if (double.IsNaN(dx) || double.IsNaN(dy))
            {
                dx = 0; dy = 0;
                Creatures[Percepteur].Vitesse = 0;
            }
            //x et y sont exprimes en % de la carte
            //dx/Scale est exprimé en case et (dx/Scale)/TerrainSubdiv est exprime en % de carte
            double x = Creatures[Percepteur].Position.Real+dx/(Scale*(double)TerrainSubdiv);
            double y = Creatures[Percepteur].Position.Imaginary + dy/(Scale * (double)TerrainSubdiv);
            if(double.IsNaN(x))
            {
                x = 0.5;
            }
            if (double.IsNaN(y))
            {
                y = 0.5;
            }
            x = Math.Min(Math.Max(float.Epsilon, x), 1.0-float.Epsilon);
            y = Math.Min(Math.Max(float.Epsilon, y), 1.0-float.Epsilon);
            if (double.IsNaN(x) || double.IsNaN(y))
            {
                x = 0.5;y = 0.5 ;
            }
            Creatures[Percepteur].Position = new Complex(x, y);
        }
        //Gestion de la population
        void PopulationUpdate()
        {
            int Morts = 0;
            //Comptage des morts
            for (int i = 0; i < EnVie.GetLength(0); i++)
            {
                if(!EnVie[i])
                {

                    Morts++;
                }
            }
            //Comptage des enfants
            int EnfantsN = Enfants.Count();
            log+=("------------> Morts : "+ Morts+"; Naissances : "+ EnfantsN+"; Solde : "+ (EnfantsN-Morts)+"\n");
            if (EnfantsN !=0 || Morts !=0)
            {
                //On update la population
                Creature[] NewPop = new Creature[EnVie.GetLength(0)-Morts+EnfantsN];
                int i;
                for (i = 0; i < Enfants.Count; i++)
                {
                    //On ajoute les enfants
                    NewPop[i] = Enfants.ElementAt(i);
                }
                for (int j = 0; j < EnVie.GetLength(0); j++)
                {
                    if (EnVie[j])
                    {
                        //On rajoute les anciens
                        NewPop[i] = Creatures[j];
                        i++;

                    }
                }
                Creatures = NewPop;
                EnVie = new bool[Creatures.GetLength(0)];
                for(int k=0;k<EnVie.GetLength(0);k++)
                {
                    //Les nouveaux individus sont en vie
                    EnVie[k] = true;
                }
                
            }
            log+=(">>>>>>>> Population Totale : " + Creatures.GetLength(0)+"<<<<<<<<"+"\n");
        }
        void AjoutIndividuAleatoire(double ProbaApparition,ref Random r, Creature MaxValues)
        {
            while (r.NextDouble()<ProbaApparition)
            {
                //On ajoute des individus aleatoires selon une certaine probabilite
                Creature C = new Creature();
                GenererCreatureAleatoire(ref C,ref r,MaxValues);
                Enfants.Add(C);
            }
        }
        //Rendu
        public void Rendu(int X,int Y, bool AffCreatures, ref Bitmap Ret)
        {
            double Intensite = 10;
            SolidBrush sb;
            Graphics G = Graphics.FromImage(Ret);
            //En pixels
            int TileWidth=(int)((double)X/(double)TerrainSubdiv);
            int TileHeight = (int)((double)Y / (double)TerrainSubdiv);
            for (int x = 0; x < TerrainSubdiv; x++)
            {
                for (int y = 0; y < TerrainSubdiv; y++)
                {
                    ElementTerrain T = Terrain[x, y];
                    //Les valeurs sont ramenees en densité sur la case
                    double Imp = T.Impraticabilité/(Scale*Scale)*Intensite;
                    double Vg = T.Vegetation/(Scale*Scale)*Intensite;
                    double Vi = T.Viande/(Scale*Scale)*Intensite;
                    int R = (int)((1.0 - Math.Exp(-Vi * 2.0 - Imp)) * 255.0);
                    int V = (int)((1.0 - Math.Exp(-Vg * 0.1 - Imp)) * 255.0);
                    int B = (int)((1.0 - Math.Exp(-Imp)) * 255.0);
                    Color clr = Color.FromArgb(R, V, B);
                    sb = new SolidBrush(clr);
                    
                    int x_ = x * TileWidth;
                    int y_ = y * TileWidth;
                    G.FillRectangle(sb,x_,y_,TileWidth,TileHeight);
                }
            }

            if (AffCreatures)
            {
                
                Font nameFont; 
                for (int i = 0; i < Creatures.GetLength(0); i++)
                {
                        int xc = (int)(Creatures[i].Position.Real * (double)X);
                        int yc = (int)(Creatures[i].Position.Imaginary * (double)Y);
                        double Rn = 2.0 * Creatures[i].Rayon/Scale; //En Cases
                        int RealRadX = Math.Max((int)(1.0 * Rn * (double)TileWidth), 0);//En pixels
                        int RealRadY = Math.Max((int)(1.0 * Rn * (double)TileWidth), 0);
                        int Rtete = (int)(RealRadX * 0.2);
                    
                    //Amange, Sestbattu,ReproducedThisTurn
                    //try
                    //{
                        Brush b = Brushes.Black;
                        if (Rn > 0 && EnVie[i] && RealRadX * 0.3 >= 1)
                        {
                            int size = (int)Math.Max((RealRadX)/5,1);
                            nameFont = new Font("Arial", size);
                            string name = Creatures[i].nom.Prenom + " " + Creatures[i].nom.NomDeFamille;
                            int C1 = Math.Min(255, Math.Max(0, (int)(Creatures[i].Couleur * 255.0)));
                            int C2 = Math.Min(255, Math.Max(0, (int)(Creatures[i].Agresivite * 255.0)));
                            int C2Bis = Math.Min(255, Math.Max(0, (int)(Sig(Creatures[i].VolonteReproductive) * 255.0)));
                            int C3 = Math.Min(255, Math.Max(0, (int)(Creatures[i].Attitude * 255.0)));

                            Complex Tete = Complex.FromPolarCoordinates(RealRadX / 2.0, Creatures[i].Vitesse.Phase);
                            //Attitude
                            sb = new SolidBrush(Color.FromArgb(C3, C3, C3));
                            G.FillEllipse(sb, xc - RealRadX / 2, yc - RealRadY / 2, RealRadX, RealRadY);
                            //Nom
                            G.DrawString(name, nameFont,b,(int)(xc+ RealRadX/3.5), yc);
                            //Repro + Agresivité
                            RealRadX = (int)((double)RealRadX * 0.7);
                            RealRadY = (int)((double)RealRadY * 0.7);
                            sb = new SolidBrush(Color.FromArgb(C2, 0, C2Bis));
                            G.FillEllipse(sb, xc - RealRadX / 2, yc - RealRadY / 2, RealRadX, RealRadY);
                            //Couleur
                            RealRadX = (int)((double)RealRadX * 0.5 / 0.7);
                            RealRadY = (int)((double)RealRadY * 0.5 / 0.7);
                            sb = new SolidBrush(Color.FromArgb(C1, C1, C1));
                            G.FillEllipse(sb, xc - RealRadX / 2, yc - RealRadY / 2, RealRadX, RealRadY);
                            //NRJ
                            RealRadX = (int)((double)RealRadX * 0.7 / 0.45);
                            RealRadY = (int)((double)RealRadY * 0.7 / 0.45);
                            Pen p = new Pen(Brushes.Green);
                            p.Width = (float)(0.02 / 0.45 * RealRadX);
                            G.DrawArc(p, xc - RealRadX / 2, yc - RealRadY / 2, RealRadX, RealRadY, 0, (float)(float)Math.Max(Math.Min(Creatures[i].Energie / 200.0, 1), 0.1) * 360);
                            RealRadX = (int)((double)RealRadX * 0.45 / 0.3);
                            RealRadY = (int)((double)RealRadY * 0.45 / 0.3);
                            p.Brush = Brushes.LightBlue;
                            G.DrawArc(p, xc - RealRadX / 2, yc - RealRadY / 2, RealRadX, RealRadY, 0, ((float)Math.Max(Creatures[i].Sante, 0.1) * 360));
                            //Tete
                            if (Creatures[i].AMange)
                            {
                                if (Creatures[i].SEstBattu)
                                {
                                    sb = new SolidBrush(Color.Yellow);
                                }
                                else
                                {
                                    sb = new SolidBrush(Color.Blue);
                                }
                            }
                            else if (Creatures[i].SEstBattu)
                            {
                                sb = new SolidBrush(Color.Red);
                            }
                            else
                            {
                                sb = new SolidBrush(Color.White);
                            }
                            G.FillEllipse(sb, xc + (int)(Tete.Real - Rtete / 2), yc + (int)(Tete.Imaginary - Rtete / 2), Rtete, Rtete);
                        }
                    //}
                   // catch(Exception e)
                    //{
                        //log+=("An error occurred: '{0}'", e);
                        //log+=(xc+" "+yc+" "+RealRadX+" "+RealRadY);
                   // }
                }
            }
            int Size =(int) (Y * 0.02);
            Font arialFont = new Font("Arial", Size);
            Brush b_ = Brushes.Black;
            G.DrawString("Temperature : "+(VE.Temperature).ToString("0.000", System.Globalization.CultureInfo.InvariantCulture), arialFont,b_,2*Size,2*Size);
            G.DrawString("Vent : " + VE.Vent.ToString("0.000", System.Globalization.CultureInfo.InvariantCulture), arialFont, b_, 2*Size, 4*Size);
            G.DrawString("Ensoleillement : " + VE.Lumiere.ToString("0.000", System.Globalization.CultureInfo.InvariantCulture), arialFont, b_, 2*Size, 6*Size);
            G.DrawString("Population : " + Creatures.GetLength(0).ToString("0", System.Globalization.CultureInfo.InvariantCulture), arialFont, b_, 2 * Size, 8 * Size);
        }
    }
    //Réseau Cerebral
    [Serializable]public class Cerveau
    {
        public int[] Structure;
        public Complex[][,] Matrices;
        double GaussianRandom(ref Random r)
        {
            double r1 = r.NextDouble();
            double r2 = r.NextDouble();
            if (r1 == 0)
            {
                r1 = 1;
            }
            double k = -Math.Log(r1);
            double rad = 2.0 * Math.PI * r2;
            return 1 + Math.Cos(rad) * k;
        }
        Complex Lissage(Complex x)
        {
            return Complex.FromPolarCoordinates(1.0-1.0/(1.0+x.Magnitude),x.Phase);
        }
        Complex[] ProduitMatricesEtLissage(Complex[,] A, Complex[] B)
        {
            Complex[] Retour = new Complex[A.GetLength(1)];
            for(int i=0;i< A.GetLength(1); i++)
            {
                Complex sum = 0;
                for(int j=0; j< A.GetLength(0); j++)
                {
                    sum += A[j, i] * B[j];
                }
                Retour[i] = Lissage(sum);
            }
            return Retour;
        }
        public void GenererAleatoire(int[] structure, ref Random r)
        {
            Structure = structure;
            Matrices = new Complex[Structure.Length - 1][,];
            for(int i =0; i<Structure.Length-1;i++)
            {
                Complex[,] Matrice = new Complex[Structure[i], Structure[i+1]];
                for(int k=0;k< Structure[i];k++)
                {
                    for (int l = 0; l < Structure[i+1]; l++)
                    {
                        Matrice[k, l] = Complex.FromPolarCoordinates(r.NextDouble(), r.NextDouble() * Math.PI * 2);
                    }
                }
                Matrices[i] = Matrice;
            }
        }
        public bool Compatible(Cerveau Partenaire)
        {
            if(Structure.Length == Partenaire.Structure.Length)
            {
                for (int i = 0; i < Structure.Length; i++)
                {
                    if (Structure[i] != Partenaire.Structure[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            else
            {
                return false;
            }
        }
        public Cerveau Moyenne(Cerveau A, Cerveau B)
        {
            Cerveau retour = new Cerveau();
            retour.Structure = A.Structure;
            retour.Matrices = new Complex[retour.Structure.Length - 1][,];
            for (int i = 0; i < retour.Matrices.Length; i++)
            {
                int x = A.Matrices[i].GetLength(0);
                int y = A.Matrices[i].GetLength(1);
                Complex[,] Matrice = new Complex[x, y];
                for (int j = 0; j < x; j++)
                {
                    for (int k = 0; k < y; k++)
                    {
                        Matrice[j, k] = (A.Matrices[i][j, k] + B.Matrices[i][j, k]) / 2;
                    }
                }
                retour.Matrices[i] = Matrice;
            }
            return retour;
        }
        public void Mutations(double TauxMutations,ref Random r)
        {
            for (int i = 0; i < Matrices.Length; i++)
            {
                int x = Matrices[i].GetLength(0);
                int y = Matrices[i].GetLength(1);
                Complex[,] Matrice = new Complex[x, y];
                for (int j = 0; j < x; j++)
                {
                    for (int k = 0; k < y; k++)
                    {
                        Matrice[j, k] = Matrices[i][j, k] +TauxMutations*Complex.FromPolarCoordinates(Math.Abs(GaussianRandom(ref r)),r.NextDouble()*Math.PI*2) ;
                    }
                }
                Matrices[i] = Matrice;
            }
        }
        public void EffectuerCalcul(Complex[] Inputs, ref Complex[] Outputs)
        {
            Complex[] Lay = Inputs;
            for (int i = 0; i < Structure.Length - 1; i++)
            {
                Lay = ProduitMatricesEtLissage(Matrices[i],Lay);
            }
            Outputs = Lay;
        }
    }
    //Generation de heightmaps
    public class Perlin
    {
        /// Perlin Noise Constructot
        public Perlin(int width, int height, ref Random r)
        {
            this.MAX_WIDTH = width;
            this.MAX_HEIGHT = height;
            GenerateNoise(ref r);
        }

        public int MAX_WIDTH = 256;
        public int MAX_HEIGHT = 256;

        /// Gets the value for a specific X and Y coordinate
        /// results in range [-1, 1] * maxHeight
        public float GenerateValues(float X, float Y, float MaxHeight,float Frequency, float Amplitude, float Persistance,int Octaves)
        {
            float FinalValue = 0.0f;
            for (int i = 0; i < Octaves; ++i)
            {
                FinalValue += GetSmoothNoise(X * Frequency, Y * Frequency) * Amplitude;
                Frequency *= 2.0f;
                Amplitude *= Persistance;
            }
            if (FinalValue < -1.0f)
            {
                FinalValue = -1.0f;
            }
            else if (FinalValue > 1.0f)
            {
                FinalValue = 1.0f;
            }
            return 0.5f*(FinalValue+1.0f) * MaxHeight;
        }

        //This function is a simple bilinear filtering function which is good (and easy) enough.        
        private float GetSmoothNoise(float X, float Y)
        {
            float FractionX = X - (int)X;
            float FractionY = Y - (int)Y;
            int X1 = ((int)X + MAX_WIDTH) % MAX_WIDTH;
            int Y1 = ((int)Y + MAX_HEIGHT) % MAX_HEIGHT;
            //for cool art deco looking images, do +1 for X2 and Y2 instead of -1...
            int X2 = ((int)X + MAX_WIDTH - 1) % MAX_WIDTH;
            int Y2 = ((int)Y + MAX_HEIGHT - 1) % MAX_HEIGHT;
            float FinalValue = 0.0f;
            FinalValue += FractionX * FractionY * Noise[X1, Y1];
            FinalValue += FractionX * (1 - FractionY) * Noise[X1, Y2];
            FinalValue += (1 - FractionX) * FractionY * Noise[X2, Y1];
            FinalValue += (1 - FractionX) * (1 - FractionY) * Noise[X2, Y2];
            return FinalValue;
        }

        float[,] Noise;
        bool NoiseInitialized = false;
        /// create a array of randoms
        private void GenerateNoise(ref Random r)
        {
            if (NoiseInitialized)                //A boolean variable in the class to make sure we only do this once
                return;
            Noise = new float[MAX_WIDTH, MAX_HEIGHT];    //Create the noise table where MAX_WIDTH and MAX_HEIGHT are set to some value>0            
            for (int x = 0; x < MAX_WIDTH; ++x)
            {
                for (int y = 0; y < MAX_HEIGHT; ++y)
                {
                    Noise[x, y] = ((float)(r.NextDouble()) - 0.5f) * 2.0f;  //Generate noise between -1 and 1
                }
            }
            NoiseInitialized = true;
        }

    }
    //Utile pour les threads
    public struct IndexInfo
    {
        //Index de depart(inclus)
        public int Start;
        //Index d'arrivee (exclus)
        public int Stop;
    }
    //Miscellianous
    [Serializable]public class Nom
    {
        public string NomDeFamille;
        public string Prenom;
        public void Generer(ref Random r)
        {
            NomDeFamille = generer(ref r);
            Prenom = generer(ref r);
        }
        public void Croiser(ref Random r, Nom P1, Nom P2, double probaChangementNomDEfamille)
        {
            Prenom = generer(ref r);
            string n1 = P1.NomDeFamille;
            string n2 = P2.NomDeFamille;
            if (r.NextDouble() < probaChangementNomDEfamille)
            {
                string nm = "";
                string[] frag = new string[3];
                frag[0] = generer(ref r);
                int r1 = r.Next(0, n1.Length - 1);
                int r2 = r.Next(0, n2.Length - 1);
                frag[1] = n1.Substring(r1, r.Next(1, n1.Length - r1));
                frag[2] = n2.Substring(r2, r.Next(1, n2.Length - r2));
                for (int i = 0; i < 3; i++)
                {
                    nm += frag[r.Next(0, 3)];
                }
                nm = nm.ToLower();
                nm = (nm[0].ToString().ToUpper()) + nm.Substring(1);
                NomDeFamille = nm;
            }
            else
            {
                if (r.NextDouble() < 0.5)
                {
                    NomDeFamille = n1;
                }
                else
                {
                    NomDeFamille = n2;
                }
            }
        }
        private string generer(ref Random r)
        {
            int lenght = r.Next(2, 10);
            string[] consonants = { "b", "c", "d", "f", "g", "h", "j", "k", "l", "m", "l", "n", "p", "q", "r", "s", "sh", "zh", "t", "v", "w", "x" ,"ph","th"};
            string[] vowels = { "a", "e", "i", "o", "u", "ae", "y", "é" };
            string Name = "";
            Name += consonants[r.Next(consonants.Length)].ToUpper();
            Name += vowels[r.Next(vowels.Length)];
            int b = 2;
            while (b < lenght)
            {
                Name += consonants[r.Next(consonants.Length)];
                b++;
                Name += vowels[r.Next(vowels.Length)];
                b++;
            }
            return Name;
        }
    }
    [Serializable]public struct ConstantesBiologiques
    {
        public double CoutManger;//En J
        public double RendementReproductionSexuee;//En %
        public double RendementClonage;//En %
        public double CoutMouvement;//En J/Effort
        public double CoutSoin;//En J/Kg/Santé
        public double CoutTemperature;// En J*Rth/°C
        public double CoutSubsistance;//En J/Kg
        public double CoutAttaque; //En J/Kg/Action
        public double MasseNRJ; //En Kg/J
        public double ExposantEffort;// En (-)
        public double IsolationGraisse; // En Rth*m
        public double DensiteGraisse; //En kg/m^3
        public double DensiteAttributs; //En Kg/m^3
        public double MasseAttributsDefenseOuAttaque; //En Kg/Action
        public double MaxViandeMangeeParKg; // en J/s/kg
        public double MaxVegetationMangeeParKg; // En Herbe/s/kg
        public double StabiliteOrganesPerception; //Sans unité
        

        public int NombreInputs;
        public int NombreOutputs;

        public ConstantesBiologiques(double coutManger, double rendementReproductionSexuee, double rendementClonage, double coutMouvement, double coutSoin, double coutTemperature, double coutSubsistance, double coutAttaque, double masseNRJ, double exposantEffort, double isolationGraisse, double densiteGraisse, double densiteAttributs, double masseAttributsDefenseOuAttaque, double maxViandeMangeeParKg, double maxVegetationMangeeParKg,double stabiliteOrganesPerception) : this()
        {
            StabiliteOrganesPerception = stabiliteOrganesPerception;
            CoutManger = coutManger;
            RendementReproductionSexuee = rendementReproductionSexuee;
            RendementClonage = rendementClonage;
            CoutMouvement = coutMouvement;
            CoutSoin = coutSoin;
            CoutTemperature = coutTemperature;
            CoutSubsistance = coutSubsistance;
            CoutAttaque = coutAttaque;
            MasseNRJ = masseNRJ;
            ExposantEffort = exposantEffort;
            IsolationGraisse = isolationGraisse;
            DensiteGraisse = densiteGraisse;
            DensiteAttributs = densiteAttributs;
            MasseAttributsDefenseOuAttaque = masseAttributsDefenseOuAttaque;
            MaxViandeMangeeParKg = maxViandeMangeeParKg;
            MaxVegetationMangeeParKg = maxVegetationMangeeParKg;

            NombreInputs = 63;
            NombreOutputs = 8;
        }
    }
    [Serializable]public struct GestionProbabilites
    {
        public Random RANDOM;
        public double FrequenceApparitionNouvelIndividu;
        public double ProbaChangementNomDeFamille;
        public float RaycastResolution;
        public double TauxMutationClonage;
        public double TauxMutationReproductionSexuee;
        public double FrequenceMaxCombat;
        public double FrequenceMaxReproduction;
        public double FrequenceMaxClonage;
        public double DiversiteLuciditeInitiale; //Sans unité, 0 = tous les sens valides pour toutes les creatures initiales, plus haut donne plus de diversite

        public GestionProbabilites(double frequenceApparitionNouvelIndividu, double probaChangementNomDeFamille, float raycastResolution, double tauxMutationClonage, double tauxMutationReproductionSexuee, double frequenceMaxCombat, double frequenceMaxReproduction, double frequenceMaxClonage, double diversiteLuciditeInitiale) : this()
        {
            DiversiteLuciditeInitiale = diversiteLuciditeInitiale;
            FrequenceApparitionNouvelIndividu = frequenceApparitionNouvelIndividu;
            ProbaChangementNomDeFamille = probaChangementNomDeFamille;
            RaycastResolution = raycastResolution;
            TauxMutationClonage = tauxMutationClonage;
            TauxMutationReproductionSexuee = tauxMutationReproductionSexuee;
            FrequenceMaxCombat = frequenceMaxCombat;
            FrequenceMaxReproduction = frequenceMaxReproduction;
            FrequenceMaxClonage = frequenceMaxClonage;
        }
    }
    [Serializable]public struct ConstantesEnvironementales
    {
        public int PremiersPeuplesNombre;
        public double Pourriture; //En 1/°C
        public double TauxDegradationViande;//En 1/s
        public double CroissanceVegetaux; //En J/s/Lumiere/Fertilite
        public double ChauffageMutuel; //En 1/m
        public double DureeJour; //En s
        public double VariabiliteVent; //En Vent
        public double VariabiliteTemp; //En °C
        public double EncombrementRelief; //En 1/m
        public double FacteurCamouflage; // En Encombrement*m
        public double SeuilMortVegetaux; // En (-)
        public double RayonInteraction; //En m
        public double ImportanceColorationVegetaux; // En m^2/herbe
        public double Gravite; // en Effort/(m/s)/kg
        public double DistanceBrouillard;

        public ConstantesEnvironementales(int premiersPeuplesNombre, double pourriture, double tauxDegradationViande, double croissanceVegetaux, double chauffageMutuel, double dureeJour, double variabiliteVent, double variabiliteTemp, double encombrementRelief, double facteurCamouflage, double seuilMortVegetaux, double rayonInteraction, double importanceColorationVegetaux, double gravite, double distanceBrouillard)
        {
            PremiersPeuplesNombre = premiersPeuplesNombre;
            Pourriture = pourriture;
            TauxDegradationViande = tauxDegradationViande;
            CroissanceVegetaux = croissanceVegetaux;
            ChauffageMutuel = chauffageMutuel;
            DureeJour = dureeJour;
            VariabiliteVent = variabiliteVent;
            VariabiliteTemp = variabiliteTemp;
            EncombrementRelief = encombrementRelief;
            FacteurCamouflage = facteurCamouflage;
            SeuilMortVegetaux = seuilMortVegetaux;
            RayonInteraction = rayonInteraction;
            ImportanceColorationVegetaux = importanceColorationVegetaux;
            Gravite = gravite;
            DistanceBrouillard = distanceBrouillard;
        }
    }
    [Serializable]public struct VariablesEnvironementales
    {
        public double Time; //En s
        public ElementTerrain[,] Terrain;
        public double Lumiere; //En Lumiere
        public double Temperature; //En °C
        public int IDlast;
        public Complex Vent; //En Vent

        public VariablesEnvironementales(double time, double lumiere, double temperature, Complex vent) : this()
        {
            Time = time;
            Lumiere = lumiere;
            Temperature = temperature;
            Vent = vent;
        }
    }
}
