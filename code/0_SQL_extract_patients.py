# Databricks notebook source
# MAGIC %md
# MAGIC #### Last updated: 03/28/2025

# COMMAND ----------

# MAGIC %md
# MAGIC # All patients in UC data warehouse

# COMMAND ----------

# MAGIC %md
# MAGIC ## Get patient demographic information

# COMMAND ----------

all_pat_info = spark.sql(
    "with gender_alt as (select distinct person_id, gender_concept_id, race_concept_id, ethnicity_concept_id, year_of_birth, location_id, concept_name as gender from omop_prod.omop_deid.person P left join omop_prod.omop_deid.concept D on P.gender_concept_id = D.concept_id), race_alt as (select person_id, gender, year_of_birth, ethnicity_concept_id, location_id, concept_name as race from gender_alt left join omop_prod.omop_deid.concept D on gender_alt.race_concept_id = D.concept_id), ethnicity_alt as (select person_id, gender, race, year_of_birth, location_id, concept_name as ethnicity from race_alt left join omop_prod.omop_deid.concept D on race_alt.ethnicity_concept_id = D.concept_id), person_clean as (select distinct person_id, gender, race, ethnicity, year_of_birth, location_source_value from ethnicity_alt inner join omop_prod.omop_deid.location on ethnicity_alt.location_id = omop_prod.omop_deid.location.location_id) select distinct person_clean.person_id, gender, race, ethnicity, year_of_birth, location_source_value, death_date from person_clean left join omop_prod.omop_deid.death on person_clean.person_id = omop_prod.omop_deid.death.person_id"
)
all_pat_info.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/all_pat_info_0328")

# COMMAND ----------

# MAGIC %md
# MAGIC Cleaned demographics patients. Extract related visit information

# COMMAND ----------

# copy files from local to dbfs
dbutils.fs.cp("file:/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/Review_Added/data/patient_data/mod/pat_clean_id_0328.csv", "s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/pat_clean_id_0328.csv")
pat_clean_id = spark.read.format("csv").option("sep",",").option("header","true").load("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/pat_clean_id_0328.csv")
pat_clean_id.createOrReplaceTempView("pat_clean_id")
pat_clean_id.show(5)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Get EHR features

# COMMAND ----------

all_pat_visit_info = spark.sql(
  "with all_pat_visit as (select distinct person_id, condition_start_date, left(concept_code, 3) as icd_3digit from omop_prod.omop_deid.condition_occurrence conditions left join omop_prod.omop_deid.concept on condition_source_concept_id = omop_prod.omop_deid.concept.concept_id where person_id in (select distinct PatientID from pat_clean_id) and vocabulary_id == 'ICD10CM') select distinct person_id, min(condition_start_date) as date_first_visit, max(condition_start_date) as date_last_visit, count(distinct condition_start_date) as n_encounter, count(distinct icd_3digit) as n_icd_3digit from all_pat_visit group by person_id"
)
all_pat_visit_info.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/all_pat_visit_0328")

# COMMAND ----------

# MAGIC %md
# MAGIC # AD patients and visits

# COMMAND ----------

# MAGIC %md
# MAGIC ## AD patient visits of all time

# COMMAND ----------

alz_pat_id = spark.sql(
  "select distinct conditions.person_id from omop_prod.omop_deid.condition_occurrence conditions left join omop_prod.omop_deid.concept on condition_source_concept_id = omop_prod.omop_deid.concept.concept_id where concept_code in ('G30', 'G30.0', 'G30.1', 'G30.8', 'G30.9')"
)
alz_pat_id.createOrReplaceTempView("alz_pat_id")
alz_pat_id.show(5)

# COMMAND ----------

# Get the number of rows
num_rows = alz_pat_id.count()
# Get the number of columns
num_cols = len(alz_pat_id.columns)
#  Print the dimensions
print("Number of rows:", num_rows) # 33899
print("Number of columns:", num_cols) # 1

# COMMAND ----------

alz_pat_visits = spark.sql(
    "with alz_pat_visit as (select distinct person_id, concept_code as icd, condition_start_date from omop_prod.omop_deid.condition_occurrence conditions left join omop_prod.omop_deid.concept on condition_source_concept_id = omop_prod.omop_deid.concept.concept_id where person_id in (select distinct person_id from alz_pat_id) and vocabulary_id == 'ICD10CM') select distinct person_id, condition_start_date, icd from (select *, row_number() over (partition by person_id, icd order by condition_start_date asc) as row_num from alz_pat_visit) subquery where row_num = 1 order by person_id asc, condition_start_date asc "
)
alz_pat_visits.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/AD_encounters_0328")

# COMMAND ----------

alz_pat_visits.show(5)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Confounding patients 

# COMMAND ----------

confound_pat_id = spark.sql(
  "select distinct conditions.person_id from omop_prod.omop_deid.condition_occurrence conditions left join omop_prod.omop_deid.concept on condition_source_concept_id = omop_prod.omop_deid.concept.concept_id where left(concept_code, 3) in ('F00', 'F01', 'F02', 'F03', 'F05', 'F06', 'F07', 'F09', 'F10', 'F18', 'F19', 'F23','F62', 'G31', 'G47', 'R40', 'R41', 'R44', 'R47', 'R48', 'R54')"
)
confound_pat_id.createOrReplaceTempView("confound_pat_id")
confound_pat_id.show(5)

# COMMAND ----------

confound_pat_id.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/confound_pat_id_0703")

# COMMAND ----------

# MAGIC %md
# MAGIC # Matched controls

# COMMAND ----------

# MAGIC %md
# MAGIC ## Direct matching

# COMMAND ----------

# copy files from local to dbfs
dbutils.fs.cp("file:/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/Review_Added/data/patient_data/mod/matched_control_id_0328.csv", "s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/matched_control_id_0328.csv")
control_pat_id = spark.read.format("csv").option("sep",",").option("header","true").load("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/matched_control_id_0328.csv")
control_pat_id.createOrReplaceTempView("control_pat_id")
control_pat_id.show(5)

# COMMAND ----------

# Get the number of rows
num_rows = control_pat_id.count()
# Get the number of columns
num_cols = len(control_pat_id.columns)
#  Print the dimensions
print("Number of rows:", num_rows) # 74472
print("Number of columns:", num_cols) # 1

# COMMAND ----------

control_pat_visits = spark.sql(
    "with control_pat_visit as (select distinct person_id, concept_code as icd, condition_start_date from omop_prod.omop_deid.condition_occurrence conditions left join omop_prod.omop_deid.concept on condition_source_concept_id = omop_prod.omop_deid.concept.concept_id where person_id in (select distinct PatientID from control_pat_id)) select distinct person_id, condition_start_date, icd from (select *, row_number() over (partition by person_id, icd order by condition_start_date asc) as row_num from control_pat_visit) subquery where row_num = 1 order by person_id asc, condition_start_date asc "
)
control_pat_visits.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/control_encounters_0328")

# COMMAND ----------

# MAGIC %md
# MAGIC ## Matching with non-confounding controls

# COMMAND ----------

# copy files from local to dbfs
dbutils.fs.cp("file:/Workspace/Users/mingzhoufu@mednet.ucla.edu/mingzhoufu@mednet.ucla.edu/AD_trajectory/data/patient_data/mod/matched_control_id_rm_0703.csv", "s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/matched_control_id_rm_0703.csv")
control_pat_id = spark.read.format("csv").option("sep",",").option("header","true").load("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/matched_control_id_rm_0703.csv")
control_pat_id.createOrReplaceTempView("control_pat_id")
control_pat_id.show(5)

# COMMAND ----------

control_pat_visits = spark.sql(
    "with control_pat_visit as (select distinct person_id, concept_code as icd, condition_start_date from omop_prod.omop_deid.condition_occurrence conditions left join omop_prod.omop_deid.concept on condition_source_concept_id = omop_prod.omop_deid.concept.concept_id where person_id in (select distinct PatientID from control_pat_id)) select distinct person_id, condition_start_date, icd from (select *, row_number() over (partition by person_id, icd order by condition_start_date asc) as row_num from control_pat_visit) subquery where row_num = 1 order by person_id asc, condition_start_date asc "
)
control_pat_visits.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("s3://uchdw-501868445017-us-west-2-prod-databricks-user-files/ucla_mingzhoufu/control_encounters_rm_0703")

# COMMAND ----------

