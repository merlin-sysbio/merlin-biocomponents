/**
 * 
 */
package pt.uminho.ceb.biosystems.merlin.biocomponents.io.readers;

import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.sbml.jsbml.xml.XMLNode;
import org.sbml.jsbml.xml.XMLTriple;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.merlin.biocomponents.container.ContainerUtils;
import pt.uminho.ceb.biosystems.merlin.database.connector.databaseAPI.CompartmentsAPI;
import pt.uminho.ceb.biosystems.merlin.database.connector.databaseAPI.ModelAPI;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.Connection;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.DatabaseAccess;
import pt.uminho.ceb.biosystems.merlin.utilities.RulesParser;
import pt.uminho.ceb.biosystems.merlin.utilities.Utilities;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.CompartmentContainer;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.MetaboliteContainer;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.model.ReactionContainer;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.CompartmentCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.GeneCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.InvalidBooleanRuleException;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.MetaboliteCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.ReactionCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.ReactionConstraintCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.ReactionTypeEnum;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.components.StoichiometryValueCI;
import pt.uminho.ceb.biosystems.mew.biocomponents.container.interfaces.IContainerBuilder;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.map.indexedhashmap.IndexedHashMap;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.pair.Pair;

/**
 * @author Oscar
 *
 */
public class ContainerBuilder implements IContainerBuilder {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(ContainerBuilder.class);

	private static final String GENE_RULE_PREFIX = "GENE_ASSOCIATION: ";
	private static final String PROTEIN_RULE_PREFIX = "PROTEIN_ASSOCIATION: ";
	private static final String PROTEIN_CLASS_PREFIX = "PROTEIN_CLASS: ";
	private static final String SUBSYSTEM_PREFIX = "SUBSYSTEM: ";
	//private static final String FOMULA_PREFIX = "FORMULA: ";
	//private static final String CHARGE_PREFIX = "CHARGE: ";

	private static final boolean concatenate = false;


	private String modelID;
	private String modelName;
	private String notes;
	private int version;
	private IndexedHashMap<String, ReactionConstraintCI> defaultEC;

	private IndexedHashMap<String, MetaboliteCI> compoundsMap;

	private IndexedHashMap<String, ReactionCI> reactionsMap;

	private IndexedHashMap<String, CompartmentCI> compartmentsMap;

	private Map<String, GeneCI> genes;

	// Extra info
	private Map<String, Map<String, String>> metabolitesExtraInfo;
	private Map<String, Map<String, String>> reactionsExtraInfo;
	private Connection connection;

	private Map<String, List<MetaboliteContainer>> reactionMetabolites;
	private Map<String, ReactionContainer> reactions;

	private Map<String, CompartmentContainer> compartments;
	private Map<String, String> compartmentID;
	private int compartmentCounter;

	private String externalCompartmentID;
	private String organismName;
	private String biomassID;

	private String biomassName;

	/**
	 * @param databaseConnector
	 * @param modelName
	 * @param isCompartmentalisedModel
	 * @param organismName
	 * @param biomassName
	 * @throws Exception
	 */
	public ContainerBuilder(DatabaseAccess databaseConnector, String modelName, boolean isCompartmentalisedModel, String organismName, String biomassName) throws Exception {

		this.connection = new Connection(databaseConnector);
		this.biomassName = biomassName;
		this.setModelID(modelID);
		this.modelName = modelName;
		this.organismName = organismName;
		this.notes = "";
		this.version = 0;
		this.compoundsMap = new IndexedHashMap<>();
		this.reactionsMap = new IndexedHashMap<>();
		this.compartmentsMap = new IndexedHashMap<>();
		this.metabolitesExtraInfo = new IndexedHashMap<>();
		this.reactionsExtraInfo	= new IndexedHashMap<>();
		this.defaultEC = new IndexedHashMap<>();
		this.reactions=new HashMap<>();
		this.reactionMetabolites=new HashMap<>();
		this.compartments=new HashMap<>();
		this.compartmentID=new HashMap<>();
		this.compartmentCounter = 1;
		//		=new HashMap<>();
		populateInformation(isCompartmentalisedModel);
	}


	/**
	 * @param isCompartmentalisedModel
	 * @throws Exception
	 */
	private void populateInformation(boolean isCompartmentalisedModel) throws Exception {

		String conditions="";

		if(isCompartmentalisedModel)
			conditions = conditions.concat(" NOT originalReaction");
		else
			conditions = conditions.concat(" originalReaction");

		Statement stmt = this.connection.createStatement();

		this.getReactions(stmt, conditions);
		this.getStoichiometry(stmt, conditions);
		this.getCompartments(stmt, isCompartmentalisedModel);
		this.buildModel();
	}


	/**
	 * @param link
	 * @param conditions
	 */
	private void getReactions(Statement stmt, String conditions) {

		try {

			stmt = connection.createStatement();

			Map<String, ArrayList<String>> result = ModelAPI.getReactions(stmt, conditions);
			ArrayList<String> list;

			for(String reactionID: result.keySet()) {

				list = result.get(reactionID);

				ReactionContainer reactionContainer = null;

				if(this.reactions.containsKey(reactionID)) {

					logger.error("same reaction {}",reactionID);
				}
				else {

					reactionContainer = new ReactionContainer(reactionID, list.get(0), list.get(1), Boolean.valueOf(list.get(2)), list.get(3), list.get(0));

					if(list.get(4)!=null)
						reactionContainer.setNotes(list.get(4));

					if(list.get(5)!= null && !list.get(5).equalsIgnoreCase("null") && !list.get(5).isEmpty())
						reactionContainer.setLowerBound(Double.parseDouble(list.get(5)));
					if(list.get(6)!= null && !list.get(6).equalsIgnoreCase("null")&& !list.get(6).isEmpty())
						reactionContainer.setUpperBound(Double.parseDouble(list.get(6)));
				}
				reactionContainer.setGeneRule(list.get(7));

				reactionContainer = this.reactions.put(reactionID, reactionContainer);
			}


			ArrayList<String[]> result2 = ModelAPI.getEnzymeHasReaction(stmt);
			String[] list2;

			for(int i=0; i<result2.size(); i++) {

				list2 = result2.get(i);	

				if(this.reactions.containsKey(list2[0]) && !list2[2].contains(".-")) {

					ReactionContainer reactionContainer = this.reactions.get(list2[0]);

					Set<String> enzymeSet = new TreeSet<String>();
					if(reactionContainer.getEnzymes()!=null) {

						enzymeSet = reactionContainer.getEnzymes();
					}

					if(list2[2]!=null) {

						enzymeSet.add(list2[2]);
						reactionContainer.setEnzymes(enzymeSet);
					}
					reactionContainer = this.reactions.put(list2[0], reactionContainer);
				}
			}


			result2 = ModelAPI.getReactionPathway(stmt);

			for(int i=0; i<result2.size(); i++) {

				list2 = result2.get(i);	

				if(this.reactions.containsKey(list2[0])) {

					ReactionContainer reactionContainer = this.reactions.get(list2[0]);

					Set<String> pathways = new TreeSet<String>();
					if(reactionContainer.getPathways()!=null) {

						pathways = reactionContainer.getPathways();
					}
					if(list2[1]!=null) {

						pathways.add(list2[2]);
						reactionContainer.setPathways(pathways);
					}
					reactionContainer = this.reactions.put(list2[0], reactionContainer);
				}
			}

			result2 = ModelAPI.getReactionGenes(stmt);

			for(int i=0; i<result2.size(); i++) {

				list2 = result2.get(i);			

				if(this.reactions.containsKey(list2[0]) && !list2[3].contains(".-")) {

					ReactionContainer reactionContainer = this.reactions.get(list2[0]);

					String locus= list2[2], geneName = null;

					if(list2[1]!=null)
						geneName = list2[1].replace(" ","").replace(",","_").replace("/","_").replace("\\","_").trim();

					this.addGeneCI(locus, geneName);
					reactionContainer.addGene(locus, geneName);
					reactionContainer = this.reactions.put(list2[0], reactionContainer);

				}

				//				System.out.println(this.geneSet);

				//						Set<String> genes = new TreeSet<String>();
				//						if(reactionContainer.getGenes()!=null) {
				//
				//							genes = reactionContainer.getGenes();
				//						}
				//
				//						if(list[1]==null || list[1].isEmpty()) {
				//
				//							genes.add(list[2].trim());
				//							reactionContainer.setGenes(genes);
				//						}
				//						else {
				//
				//							genes.add(list[1].replace(" ","").replace(",","_").replace("/","_").replace("\\","_")
				//									//.replace("-","_")
				//									.trim()+"_"+list[2].trim());
				//							reactionContainer.setGenes(genes);
				//						}
			}

			stmt.close();
		}
		catch (SQLException e) {

			e.printStackTrace();
		}
	}

	private void addGeneCI(String locus, String geneName) {

		if(this.genes==null)
			this.genes = new IndexedHashMap<>();

		if(!this.genes.containsKey(locus))
			this.genes.put(locus, new GeneCI(locus, geneName));
	}


	/**
	 * @param link
	 * @param conditions
	 */
	private void getStoichiometry(Statement stmt, String conditions) {

		try {

			ArrayList<String[]> result2 = ModelAPI.getStoichiometryInfo(stmt, conditions);
			String[] list2;

			for(int i=0; i<result2.size(); i++) {

				list2 = result2.get(i);			

				if(this.reactions.containsKey(list2[1])) {

					if(!list2[4].contains("m") && !list2[4].contains("n")) {

						List<MetaboliteContainer> metabolitesContainer = new ArrayList<MetaboliteContainer>();

						if(this.reactionMetabolites.containsKey(list2[1])) {

							metabolitesContainer = this.reactionMetabolites.get(list2[1]);
						}

						MetaboliteContainer metabolite = new MetaboliteContainer(Integer.parseInt(list2[2]), list2[6], list2[7], list2[4], list2[5], list2[3]);
						metabolite.setEntryID(list2[8]);
						metabolitesContainer.add(metabolite);

						this.reactionMetabolites.put(list2[1], metabolitesContainer);
					}
					else {

						this.reactionMetabolites.remove(list2[1]);
						this.reactions.remove(list2[1]);
					}
				}
			}

		} catch (SQLException e) {
			e.printStackTrace();
		}
	}


	/**
	 * @param stmt
	 * @param isCompartmentalized
	 */
	private void getCompartments(Statement stmt, boolean isCompartmentalized) {

		try {

			stmt = connection.createStatement();

			Map<String, ArrayList<String>> result = CompartmentsAPI.getCompartmentsInfo(stmt);
			ArrayList<String> list;

			for(String idCompartment : result.keySet()) {

				list = result.get(idCompartment);

				CompartmentContainer compartmentContainer = new CompartmentContainer(idCompartment, list.get(0), list.get(1));
				this.compartments.put(idCompartment, compartmentContainer);

				if((list.get(0).equalsIgnoreCase("extracellular") && isCompartmentalized) || (list.get(0).equalsIgnoreCase("outside") && !isCompartmentalized))
					this.externalCompartmentID=this.getCompartmentID(idCompartment);
			}

		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param compartment
	 * @return
	 */
	private String getCompartmentID(String compartment) {

		if(!this.compartmentID.containsKey(compartment)) {

			String id = ContainerBuilder.buildID("C_", this.compartmentCounter);

			CompartmentCI c = new CompartmentCI(id, this.compartments.get(compartment).getName(), this.externalCompartmentID);
			this.compartmentsMap.put(id, c);
			this.compartmentID.put(compartment,id);
			this.compartmentCounter++ ;

		}
		return this.compartmentID.get(compartment);
	}

	/**
	 * @param type
	 * @param counter
	 * @return
	 */
	private static String buildID(String type, int counter) {

		if(counter<10 ){return type.concat("0000"+counter);}
		if(counter<100) {return type.concat("000"+counter);}
		if(counter<1000) {return type.concat("00"+counter);}
		if(counter<10000) {return type.concat("0"+counter);}

		return type.concat("_"+counter);
	}


	/**
	 * @throws InvalidBooleanRuleException 
	 * @throws SQLException 
	 * @throws NumberFormatException 
	 * 
	 */
	private void buildModel() throws InvalidBooleanRuleException, NumberFormatException, SQLException {

		int reactionsCounter = 1 ;
		int metabolitesCounter = 1;

		Map<String,String> compoundCompartmentID = new TreeMap<String, String>();

		Statement stmt = connection.createStatement();

		for(String reaction_id : this.reactions.keySet()) {

			ReactionContainer reaction = this.reactions.get(reaction_id);
			String name = reaction.getName()+"__("+reaction.getEquation().replace(" ", "")+")"+"__"+reaction_id;

			String rid = ContainerBuilder.buildID("R_", reactionsCounter)/*+"_"+reaction.getName().replace(" ", "_").replace("\t", "_").replace("-", "_")*/;
			//			String rid = reaction_id;
			reactionsCounter++ ;

			if(this.biomassName!=null && this.biomassName.equalsIgnoreCase(reaction.getName()))
				this.biomassID = rid;

			Map<String, StoichiometryValueCI> reactants = new HashMap<String, StoichiometryValueCI>();
			Map<String, StoichiometryValueCI> products = new HashMap<String, StoichiometryValueCI>();
			if(this.reactionMetabolites.get(reaction_id) != null) {

				MetaboliteCI mci = null;
				for(MetaboliteContainer metabolite : this.reactionMetabolites.get(reaction_id)) {

					String metabolite_surrogate = metabolite.getMetaboliteID()+"".concat("_").concat(metabolite.getCompartment_name());
					String mid;
					if(compoundCompartmentID.containsKey(metabolite_surrogate)) {

						mid = compoundCompartmentID.get(metabolite_surrogate);
						mci = this.compoundsMap.get(mid);
					}
					else {

						mid = ContainerBuilder.buildID("M_", metabolitesCounter);
						//						mid = Integer.toString(metabolite.getMetaboliteID());
						metabolitesCounter++ ;
						compoundCompartmentID.put(metabolite_surrogate,mid);

						if(this.compoundsMap.containsKey(mid)) {

							mci = this.compoundsMap.get(mid);
						}
						else  {

							String containerMetaboliteName="";

							if(metabolite.getEntryID() != null) {

								containerMetaboliteName=containerMetaboliteName.concat(metabolite.getEntryID());
							}

							if(metabolite.getName() != null) {

								if(containerMetaboliteName!=null && !containerMetaboliteName.isEmpty()) {

									containerMetaboliteName=containerMetaboliteName.concat("_");
								}
								containerMetaboliteName=containerMetaboliteName.concat(metabolite.getName());
							}

							if(metabolite.getFormula() != null) {

								if(containerMetaboliteName!=null && !containerMetaboliteName.isEmpty()) {

									containerMetaboliteName=containerMetaboliteName.concat("_");
								}
								containerMetaboliteName=containerMetaboliteName.concat(metabolite.getFormula());
							}

							mci = new MetaboliteCI(mid, containerMetaboliteName);

							Map<String, String> value = new HashMap<>();
							value.put("KEGG_CPD", metabolite.getEntryID());
							value.put("MERLIN_ID", metabolite.getMetaboliteID()+"");

							this.metabolitesExtraInfo.put(mid, value);

							if(metabolite.getFormula() != null)
								mci.setFormula(metabolite.getFormula());
							this.compoundsMap.put(mid, mci);
						}
					}

					mci.addReaction(reaction_id);

					String compartmentId = this.getCompartmentID(metabolite.getCompartment_name());

					if(compartmentId==null)
						logger.error("null compartment {} for metabolite {} ", metabolite,this.compartmentID);

					double value = Double.valueOf(metabolite.getStoichiometric_coefficient());

					StoichiometryValueCI s = new StoichiometryValueCI(mid, Math.abs( value), compartmentId);
					this.compartmentsMap.get(compartmentId).addMetaboliteInCompartment(mid);
					if(value>0)
						products.put(mid, s);
					else
						reactants.put(mid, s);
				}
			}
			else {

				logger.error("null reaction metabolites for {} {}",reaction_id ,reaction.getEntryID());
			}

			double upper_bound = 999999;
			double lower_bound = 0;

			if(reaction.getLowerBound()!= null)
				lower_bound = reaction.getLowerBound();
			else 
				if(reaction.isReversible())
					lower_bound = -999999;

			if(reaction.getUpperBound()!= null)
				upper_bound = reaction.getUpperBound();

			this.defaultEC.put(rid, new ReactionConstraintCI(lower_bound, upper_bound));

			ReactionCI r = new ReactionCI(rid, name, reaction.isReversible(), reactants, products);
			if(this.reactions.get(reaction_id).getEnzymes()!=null)
				r.setEc_number(this.reactions.get(reaction_id).getEnzymes().toString());

			if(ContainerUtils.isReactionDrain(r))
				r.setType(ReactionTypeEnum.Drain);

			//ADD PATHWAYS
			if(reaction.getPathways()!=null && !reaction.getPathways().isEmpty())
				r.setSubsystems(reaction.getPathways());

			//ADD GENES
			if(reaction.getGenes()!= null)
				for(Pair<String,String> gene : reaction.getGenes())			
					r.addGene(gene.getA());

			String geneRule = null;

			if(reaction.getGeneRule()==null || reaction.getGeneRule().isEmpty())
				geneRule = RulesParser.processReactionGenes(reaction.getGenes(), concatenate);

			else {

				List<List<Pair<String,String>>> rule = ModelAPI.parseBooleanRule(reaction.getGeneRule(), stmt);
				geneRule = Utilities.parseRuleListToString(rule);
			}

//			String gRule = ContainerBuilder.getSimpleGeneRule(ContainerBuilder.processReactionGenes(reaction.getGenes()), reaction.getGeneRule());
			//			String gRule = MerlinDBReader.getSimpleGeneRule_SBML2(MerlinDBReader.processReactionGenes(reaction.getGenes()), MerlinDBReader.processReactionNotes(reaction.getNotes()));
			
			
			if(geneRule!= null && !geneRule.isEmpty())
				r.setGeneRule(geneRule);

			//				else
			//					System.out.println(reaction_id+"\t"+r.getName().split("__")[0]);


			//			}

			//			String genesNotes = MerlinDBReader.processReactionGenes(reaction.getGenes());
			//			String proteinsNotes = MerlinDBReader.processReactionProteins(reaction.getEnzymes());
			//			String pathwaysNotes = MerlinDBReader.processReactionPathways(reaction.getPathways());
			//			Set<String> notesList = MerlinDBReader.processReactionNotes(notes);

			//			Set <XMLNode> g = MerlinDBReader.getGeneRules(genesNotes, notesList);
			//			
			//			if(g!=null)
			//				r.setGeneRule(g.toString());
			//			
			//			Set <XMLNode> pr = MerlinDBReader.getProteinRules(proteinsNotes, notesList);
			//			if(pr!=null)
			//				r.setProteinRule(pr.toString().replace(".", "").replace("#", ""));
			//
			//			Set <XMLNode> pa = MerlinDBReader.getPathwaysRules(pathwaysNotes, notesList);
			//			if(pa!=null)
			//				r.setSubsystem(pa.toString());

			this.reactionsMap.put(r.getId(), r);

			Map<String, String> value = new HashMap<>();
			value.put("MERLIN_ID", reaction.getEntryID());

			this.reactionsExtraInfo.put(rid, value);
		}

		stmt.close();
	}


	/**
	 * @param proteins
	 * @return
	 */
	public static String processReactionProteins(Set<String> proteins) {

		String proteinData = "";

		if(proteins!=null) {

			for(String prot:proteins) {

				if(prot!=null) {

					proteinData=proteinData.concat(prot).concat(" or ");
				}
			}
			proteinData=proteinData.substring(0, proteinData.lastIndexOf(" or "));
		}
		return proteinData;
	}

	/**
	 * @param pathways
	 * @return
	 */
	public static String processReactionPathways(Set<String> pathways) {

		String pathwaysData = "";

		if(pathways!=null) {

			for(String path:pathways) {

				if(path!=null) {

					pathwaysData=pathwaysData.concat(path).concat(" and ");
				}
			}
			pathwaysData=pathwaysData.substring(0, pathwaysData.lastIndexOf(" and "));
		}
		return pathwaysData;
	}

	/**
	 * @param enzymes
	 * @return
	 */
	public static String processReactionProteinClass(Set<String> enzymes) {

		String enzymesData = "";

		if(enzymes!=null) {

			for(String enzyme:enzymes) {

				if(enzyme!=null)
					enzymesData=enzymesData.concat(enzyme).concat(" and ");
			}
			enzymesData=enzymesData.substring(0, enzymesData.lastIndexOf(" and "));
		}
		return enzymesData;	
	}


	/**
	 * @param notes
	 * @return
	 */
	public static Set<String> processReactionNotes(String notes) {

		Set<String> reactionNotes = new HashSet<>();
		//<html:p>GENE_ASSOCIATION: </html:p>
		//<html:p>PROTEIN_ASSOCIATION: "</html:p>
		//<html:p>SUBSYSTEM: S_Tyrosine__Tryptophan__and_Phenylalanine_Metabolism</html:p>
		//<html:p>PROTEIN_CLASS: 1.13.11.27</html:p>
		//<html:p>GENE_ASSOCIATION: ( YOL096C  and  YDR204W  and  YML110C  and  YGR255C  and  YOR125C  and  YGL119W  and  YLR201C )</html:p>
		//<html:p>PROTEIN_ASSOCIATION: ( Coq3-m and Coq4-m and Coq5-m and Coq6-m and Coq7-m and Coq8-m and Coq9-m )"</html:p>
		//<html:p>SUBSYSTEM: S_Quinone_Biosynthesis</html:p>
		//<html:p>PROTEIN_CLASS: </html:p>


		if(notes!=null && !notes.trim().isEmpty()) {

			for(String note : notes.split(" \\| "))
				reactionNotes.add(note.replace(",","_").replace(" ","_").replace(":_",": ").replace("_AND_"," and ").replace("_OR_"," or ").trim());	 //.replace(")","").replace("(","_")
		}
		return reactionNotes; 
	}

	/**
	 * @param genesNotes
	 * @param notesList
	 * @param addAllNotes
	 * @param gene_rule
	 * @return
	 */
	public static Set<XMLNode> getGeneRules(String merlinGenes, Set<String> notesList, boolean addAllNotes, String gene_rule) {

		Set<XMLNode> xmlNodeSet = new HashSet<>();
		XMLNode node = new XMLNode(new XMLTriple("p", "", "html"));
		boolean noGene = true, go = true;
		String genesNotes = null;

		if(gene_rule!=null && !gene_rule.isEmpty()) {

			noGene = false;
			genesNotes = gene_rule;
		}

		for (String noteReaction : notesList) {

			if(noteReaction.contains(GENE_RULE_PREFIX) && genesNotes==null) {

				go = false;
				noGene = false;

				if(noteReaction.contains("no_gene")) {

					if(addAllNotes) {

						node = new XMLNode(new XMLTriple("p", "", "html"));
						node.addChild(new XMLNode("reaction annotation: no gene"));
						xmlNodeSet.add(node);
					}
				}
				else {

					node = new XMLNode(new XMLTriple("p", "", "html"));

					//					String out = GENE_RULE_PREFIX;
					//					
					//					String note = noteReaction.split(":")[1];
					//					
					//					String [] rules = note.split(" or ");
					//					
					//					String auxGeneOR = " or ";
					//					
					//					for(int j = 0 ; j<rules.length; j++) {
					//						
					//						String rule = rules[j];
					//						
					//						String [] genes = rule.split(" and ");
					//						
					//						String auxGeneAND = " and ";
					//						
					//						for(int i = 0 ; i<genes.length; i++) {
					//							
					//							String gene = genes[i];
					//							
					//							String[] gene_locus = gene.split("_");
					//							
					//							if(genes.length-1 == i)
					//								auxGeneAND = "";
					//							
					//							out = out.concat(gene_locus[gene_locus.length-1]).concat(auxGeneAND);
					//						}
					//						
					//						if(j<rules.length-1)
					//							out = out.concat(auxGeneOR);
					//					}
					//					
					//					node.addChild(new XMLNode(out.replace(" OR ", " or ").replace(" AND ", " and ")));

					node.addChild(new XMLNode(noteReaction.replace(" OR ", " or ").replace(" AND ", " and ")));
					xmlNodeSet.add(node);
				}
			}

			if(addAllNotes && go) {

				node = new XMLNode(new XMLTriple("p", "", "html"));
				node.addChild(new XMLNode(noteReaction));
				xmlNodeSet.add(node);
			}
		}


		if(noGene) {

			genesNotes = merlinGenes;
		}
		else {

			if(merlinGenes!=null && !merlinGenes.trim().isEmpty() && addAllNotes) {

				node = new XMLNode(new XMLTriple("p", "", "html"));
				node.addChild(new XMLNode("merlin_gene_data: "+merlinGenes));
				xmlNodeSet.add(node);
			}
		}

		if(go) {

			String ret="";
			if(genesNotes!=null && !genesNotes.trim().isEmpty())
				ret += genesNotes;
			node = new XMLNode(new XMLTriple("p", "", "html"));
			node.addChild(new XMLNode(GENE_RULE_PREFIX+ret));
			xmlNodeSet.add(node);
		}

		return xmlNodeSet;		
	}

	/**
	 * This method returns the gene rule as a String.
	 * 
	 * @param genesNotes
	 * @param notesList
	 * @param addAllNotes 
	 * @return
	 */
	public static String getSimpleGeneRule_SBML2(String genesNotes, Set<String> notesList) {

		String ret= null;
		boolean noGene = true;

		for (String noteReaction : notesList) {

			if(noteReaction.contains(GENE_RULE_PREFIX)) {

				if(!noteReaction.contains("no_gene")) 

					noGene = false;

				ret = noteReaction.replace(" OR ", " or ").replace(" AND ", " and ").replace(GENE_RULE_PREFIX, "").trim();
			}
		}

		if(noGene) {

			if(genesNotes!=null && !genesNotes.trim().isEmpty())
				ret = genesNotes;
		}

		return ret;		
	}

	/**
	 * This method returns the gene rule as a String.
	 * 
	 * @param genesNotes
	 * @param boolean_rule 
	 * @return
	 */
	public static String getSimpleGeneRule(String genesNotes, String boolean_rule) {

		String ret= null;

		if (boolean_rule != null && !boolean_rule.isEmpty())
			ret = boolean_rule;
		else if(genesNotes!=null && !genesNotes.trim().isEmpty())
			ret = genesNotes;

		return ret;		
	}

	/**
	 * @param proteinNotes
	 * @param notesList
	 * @param addAllNotes 
	 * @return
	 */
	public static Set<XMLNode> getProteinRules(String proteinNotes, Set<String> notesList, boolean addAllNotes) {

		Set<XMLNode> xmlNodeSet = new HashSet<>(); 
		XMLNode node = new XMLNode(new XMLTriple("p", "", "html"));
		boolean noProtein = true;

		for (String noteReaction : notesList) {

			if(noteReaction.contains(PROTEIN_RULE_PREFIX)) {

				noProtein = false;

				node = new XMLNode(new XMLTriple("p", "", "html"));
				node.addChild(new XMLNode(noteReaction.replace(" OR ", " or ").replace(" AND ", " and ")));
				xmlNodeSet.add(node);

				if(proteinNotes!=null && !proteinNotes.trim().isEmpty() && addAllNotes) {

					node = new XMLNode(new XMLTriple("p", "", "html"));
					node.addChild(new XMLNode("merlin_protein_data: "+proteinNotes));
					xmlNodeSet.add(node);
				}
			}
			//			else {
			//
			//				node.addChild(new XMLNode(noteReaction));
			//			}
		}

		if(noProtein) {

			String ret = "";
			if(proteinNotes!=null && !proteinNotes.trim().isEmpty())
				ret += proteinNotes;

			node = new XMLNode(new XMLTriple("p", "", "html"));
			node.addChild(new XMLNode(PROTEIN_RULE_PREFIX+ret));
			xmlNodeSet.add(node);
		}

		return xmlNodeSet;	
	}

	/**
	 * @param pathwayNotes
	 * @param notesList
	 * @param addAllNotes 
	 * @return
	 */
	public static Set<XMLNode> getPathwaysRules(String pathwayNotes, Set<String> notesList, boolean addAllNotes) {

		Set<XMLNode> xmlNodeSet = new HashSet<>();
		XMLNode node = new XMLNode(new XMLTriple("p", "", "html"));
		boolean noPathway = true;

		for (String note : notesList) {

			if(note.contains(SUBSYSTEM_PREFIX)) {

				noPathway = false;

				node = new XMLNode(new XMLTriple("p", "", "html"));
				node.addChild(new XMLNode(note.replace(" OR ", " or ").replace(" AND ", " and ")));
				xmlNodeSet.add(node);

				if(pathwayNotes!=null && !pathwayNotes.trim().isEmpty() && addAllNotes) {

					node = new XMLNode(new XMLTriple("p", "", "html"));
					node.addChild(new XMLNode("merlin_subsystem_data: "+pathwayNotes));
					xmlNodeSet.add(node);
				}
			}
		}

		if(noPathway) {

			String ret = "";
			if(pathwayNotes!=null && !pathwayNotes.trim().isEmpty())
				ret += pathwayNotes;

			node = new XMLNode(new XMLTriple("p", "", "html"));
			node.addChild(new XMLNode(SUBSYSTEM_PREFIX+ret));
			xmlNodeSet.add(node);
		}

		return xmlNodeSet;	
	}

	/**
	 * @param enzymesNotes
	 * @param notesList
	 * @param addAllNotes 
	 * @return
	 */
	public static Set<XMLNode> getProteinClassRules(String enzymesNotes, Set<String> notesList, boolean addAllNotes) {

		Set<XMLNode> xmlNodeSet = new HashSet<>();
		XMLNode node = new XMLNode(new XMLTriple("p", "", "html"));
		boolean noEnzyme = true;

		for (String note : notesList) {

			if(note.contains(PROTEIN_CLASS_PREFIX)) {

				noEnzyme = false;
				node = new XMLNode(new XMLTriple("p", "", "html"));
				node.addChild(new XMLNode(note.replace(" OR ", " or ").replace(" AND ", " and ")));
				xmlNodeSet.add(node);

				if(enzymesNotes!=null && !enzymesNotes.trim().isEmpty() && addAllNotes) {

					node = new XMLNode(new XMLTriple("p", "", "html"));
					node.addChild(new XMLNode("merlin_protein_class_data: "+enzymesNotes));
					xmlNodeSet.add(node);
				}
			}
		}
		if(noEnzyme) {

			String ret = "";
			if(enzymesNotes!=null && !enzymesNotes.trim().isEmpty())
				ret += enzymesNotes;

			node = new XMLNode(new XMLTriple("p", "", "html"));
			node.addChild(new XMLNode(PROTEIN_CLASS_PREFIX+ret));
			xmlNodeSet.add(node);
		}

		return xmlNodeSet;	
	}

	/**
	 * @return The model name
	 */
	@Override
	public String getModelName() {
		return modelName;
	}

	/**
	 * @return The organism name
	 */
	@Override
	public String getOrganismName() {
		return organismName;
	}

	/**
	 * @return The notes
	 */
	@Override
	public String getNotes() {
		return notes;
	}

	/**
	 * @return The version
	 */
	@Override
	public Integer getVersion() {
		return version;
	}

	/**
	 * @return A map with the reactions
	 */
	@Override
	public Map<String, ReactionCI> getReactions() {
		return reactionsMap;
	}

	/**
	 * @return A map with the metabolites
	 */
	@Override
	public Map<String, MetaboliteCI> getMetabolites() {
		return compoundsMap;
	}

	/**
	 * @return A map with the genes
	 */
	@Override
	public Map<String, GeneCI> getGenes() {
		return genes;
	}

	/**
	 * @return A map with the compartments
	 */
	@Override
	public HashMap<String, CompartmentCI> getCompartments() {
		return compartmentsMap;
	}

	/**
	 * @return The biomassID
	 */
	@Override
	public String getBiomassId() {
		return biomassID;
	}

	/**
	 * @return A map with the Environmental Conditions
	 */
	@Override
	public Map<String, ReactionConstraintCI> getDefaultEC() {
		return defaultEC;
	}

	/**
	 * @return The external compartment ID
	 */
	@Override
	public String getExternalCompartmentId() {
		return this.externalCompartmentID;
	}

	/**
	 * @return A map with the metabolites extra info
	 */
	@Override
	public Map<String, Map<String, String>> getMetabolitesExtraInfo() {
		return metabolitesExtraInfo;
	}

	/**
	 * @return A map with the reactions extra info
	 */
	@Override
	public Map<String, Map<String, String>> getReactionsExtraInfo() {
		return reactionsExtraInfo;
	}

	public String getModelID() {
		return modelID;
	}

	public void setModelID(String modelID) {
		this.modelID = modelID;
	}
}
