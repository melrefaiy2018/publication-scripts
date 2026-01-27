# Zenodo Upload Guide for IsiA Paper Repository

This guide explains how to upload the IsiA_paper repository to Zenodo and how each documentation file will be used.

---

## Overview

Zenodo is a free, open-access repository for research data and software. This guide covers two upload methods and explains what happens to each file during the process.

---

## Upload Methods

### Method 1: GitHub Integration (Recommended)

This method automatically syncs your GitHub repository with Zenodo.

#### Step 1: Enable GitHub-Zenodo Integration

1. **Go to Zenodo:** https://zenodo.org/
2. **Log in** using your GitHub account (or link your GitHub account if using email login)
3. **Go to your profile** → Click your username in top right
4. **Click "GitHub"** in the dropdown menu
5. **Enable the repository:**
   - Find "Papers" or "publication-scripts" repository in the list
   - Toggle the switch to **ON** (green)
6. **Sync repositories** if you don't see your repo (click the Sync button)

#### Step 2: Create a GitHub Release

1. **Go to your GitHub repository:**
   ```
   https://github.com/melrefaiy2018/publication-scripts
   ```

2. **Create a new release:**
   - Click "Releases" (right sidebar)
   - Click "Create a new release"
   - **Tag version:** `v1.0.0` (or `2025.01` for date-based versioning)
   - **Release title:** "Data and Code for IsiA Quenching Study"
   - **Description:** Brief summary (can copy from .zenodo.json description)
   - Click "Publish release"

3. **Zenodo automatically creates a record:**
   - Within minutes, Zenodo will capture the release
   - It creates a permanent DOI
   - It archives the entire repository

#### Step 3: Verify and Publish on Zenodo

1. **Go to Zenodo:** https://zenodo.org/account/settings/github/
2. **Find your repository** in the list (should show "Enabled")
3. **Click the repository name** to see the Zenodo record
4. **Verify the metadata:**
   - ✓ Title from .zenodo.json should be auto-filled
   - ✓ Authors from .zenodo.json should be listed
   - ✓ Description should be populated
   - ✓ Keywords should be present
   - ✓ License should be MIT

5. **Edit if needed:**
   - Click "Edit" to modify any fields
   - Add additional metadata (funding, references, etc.)
   - Add communities (optional - for discovery)

6. **Publish:**
   - Click "Publish" to make it public
   - DOI is now permanent and citable

---

### Method 2: Direct Upload (Alternative)

If you prefer not to use GitHub integration:

#### Step 1: Prepare Archive

```bash
cd /Users/mohamed/Documents/Research/Projects/GitHub_repo/Papers/2025/

# Create a clean archive (without .git directory)
tar -czf IsiA_paper_zenodo.tar.gz \
    --exclude='IsiA_paper/.git' \
    --exclude='IsiA_paper/__pycache__' \
    --exclude='IsiA_paper/**/__pycache__' \
    --exclude='IsiA_paper/*.pyc' \
    --exclude='IsiA_paper/**/*.pyc' \
    --exclude='IsiA_paper/.DS_Store' \
    --exclude='IsiA_paper/**/.DS_Store' \
    IsiA_paper/

# Or create a ZIP archive
zip -r IsiA_paper_zenodo.zip IsiA_paper/ \
    -x "IsiA_paper/.git/*" \
    -x "IsiA_paper/__pycache__/*" \
    -x "IsiA_paper/**/__pycache__/*" \
    -x "IsiA_paper/*.pyc" \
    -x "IsiA_paper/.DS_Store"
```

#### Step 2: Upload to Zenodo

1. **Go to Zenodo:** https://zenodo.org/
2. **Log in** to your account
3. **Click "New Upload"** (top right or from your dashboard)
4. **Upload your archive:**
   - Drag and drop `IsiA_paper_zenodo.tar.gz` or `IsiA_paper_zenodo.zip`
   - Wait for upload to complete
5. **Check if .zenodo.json was read:**
   - If present in the archive, Zenodo should auto-fill metadata
   - If not, manually fill in the form

#### Step 3: Fill in Metadata Form

If .zenodo.json wasn't auto-detected, manually enter:

- **Upload type:** Dataset
- **Title:** "Data and Code for: Quenching of the Photosynthetic Antenna IsiA is Facilitated by its Red-Emitting States"
- **Authors:** (Copy from .zenodo.json)
  - Mohamed A. A. Elrefaiy (MIT)
  - Dvir Harris
  - Hila Toporik
  - Christopher J. Gisriel
  - Yuval Mazor
  - Doran I. G. B. Raccah
  - Gabriela S. Schlau-Cohen (MIT)
- **Description:** (Copy from .zenodo.json)
- **License:** MIT License
- **Keywords:** photosynthesis, spectroscopy, protein dynamics, hamiltonian, fluorescence, computational chemistry

#### Step 4: Additional Metadata (Optional)

- **Related identifiers:**
  - Add manuscript DOI when published (relation: "is supplemented by this upload")
- **Contributors:** (if different from authors)
- **Funding:** Add grant numbers (from ACKNOWLEDGMENTS.md)
- **Communities:** Search for relevant communities (e.g., "Photosynthesis", "Computational Chemistry")

#### Step 5: Publish

- Click "Publish" to make it public
- DOI is now permanent

---

## What Happens to Each File on Zenodo?

### Files Zenodo Uses Automatically

| File | How Zenodo Uses It | Visibility |
|------|-------------------|------------|
| `.zenodo.json` | **Auto-fills upload form** with metadata | Not displayed (metadata only) |
| `CITATION.cff` | **Backup metadata source** (if .zenodo.json missing) | Can be displayed as a file |
| `README.md` | **Displayed on landing page** (if Zenodo detects it) | ✓ Visible |
| `LICENSE` | **Shown in license section** | ✓ Visible as file |

### Files for User Documentation (Included in Archive)

| File | Purpose on Zenodo | Visibility |
|------|------------------|------------|
| `AUTHORS.md` | Available for download, human-readable attribution | ✓ Downloadable |
| `CHANGELOG.md` | Shows version history to users | ✓ Downloadable |
| `QUICKSTART.md` | Helps users get started quickly | ✓ Downloadable |
| `ACKNOWLEDGMENTS.md` | Funding and credit information | ✓ Downloadable |
| All other files | Part of the research archive | ✓ Downloadable |

---

## After Publishing on Zenodo

### Step 1: Get Your DOI

After publishing, Zenodo assigns a permanent DOI, for example:
```
https://doi.org/10.5281/zenodo.1234567
```

### Step 2: Update Repository with DOI

1. **Update README.md badge:**

Find this line in README.md (currently has placeholder):
```markdown
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXXX-blue)](https://doi.org/10.5281/zenodo.XXXXXX)
```

Replace with actual DOI:
```markdown
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.1234567-blue)](https://doi.org/10.5281/zenodo.1234567)
```

2. **Commit the change:**
```bash
cd /Users/mohamed/Documents/Research/Projects/GitHub_repo/Papers/2025/IsiA_paper
git add README.md
git commit -m "Update README with Zenodo DOI badge

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
git push
```

### Step 3: Reference in Your Paper

Add to your manuscript's "Data Availability" section:
```
The data and analysis code for this study are openly available
on Zenodo at https://doi.org/10.5281/zenodo.1234567
```

---

## Versioning on Zenodo

### Creating New Versions

If you need to update your data/code after initial publication:

1. **Create a new GitHub release:**
   - Tag: `v1.1.0` (increment version)
   - Zenodo automatically creates a new version
   - **New DOI for new version** + **Concept DOI for all versions**

2. **Or upload new version directly on Zenodo:**
   - Go to your published record
   - Click "New version"
   - Upload updated files
   - Zenodo links versions together

### DOI Structure
- **Concept DOI:** Points to latest version (e.g., 10.5281/zenodo.1234567)
- **Version DOIs:** Points to specific versions (e.g., 10.5281/zenodo.1234568, 10.5281/zenodo.1234569)

**Recommendation:** Use concept DOI in publications (always points to latest)

---

## Troubleshooting

### .zenodo.json Not Working

**Symptoms:** Zenodo doesn't auto-fill metadata

**Solutions:**
1. Check JSON syntax: Use a JSON validator (https://jsonlint.com/)
2. Ensure .zenodo.json is in the **root directory** of your repository
3. Try direct upload method instead of GitHub integration
4. Manually fill in the form using .zenodo.json as reference

### GitHub Integration Not Showing Repository

**Solutions:**
1. Click "Sync now" button on Zenodo's GitHub page
2. Make sure repository is public (Zenodo can't access private repos without special permissions)
3. Check that Zenodo app has permission to access your GitHub account
4. Try disconnecting and reconnecting GitHub integration

### Upload Too Large

**If repository is very large (>50 GB):**
1. Consider splitting into multiple datasets
2. Use Zenodo's file-by-file upload
3. Consider excluding large intermediate files (keep only essential data)
4. Contact Zenodo support for large datasets (they can accommodate >50 GB)

### Missing Files After Upload

**Check .gitignore:**
- Some files might be excluded by .gitignore
- For Zenodo upload, you may want to include some ignored files (like .npy/.npz data)
- Use tar/zip commands above that explicitly include needed files

---

## Best Practices

### 1. Version Your Release Properly
- Use semantic versioning: v1.0.0 (major.minor.patch)
- Or date-based: 2025.01 (year.month)
- Tag releases in chronological order

### 2. Write a Good Description
Your .zenodo.json description is the first thing people see:
- Clearly state what the dataset contains
- Mention key results or use cases
- Reference the associated paper

### 3. Choose Relevant Communities
Browse Zenodo communities and join 1-3 relevant ones:
- Increases discoverability
- Connects your work to similar research
- Examples: "Photosynthesis", "Structural Biology", "Computational Chemistry"

### 4. Add Funding Information
- Include grant numbers in Zenodo metadata
- Helps track research impact
- Required by some funding agencies

### 5. Link to Your Paper
Once your manuscript is published:
- Update .zenodo.json "related_identifiers" with paper DOI
- Create new version on Zenodo if needed
- Update GitHub README with paper link

---

## Example: Complete Upload Workflow

### Timeline

**Before Submission:**
1. ✓ All documentation files created (AUTHORS.md, CITATION.cff, etc.)
2. ✓ .zenodo.json prepared with placeholder manuscript DOI
3. ✓ Repository tested and verified
4. ✓ .gitignore properly configured

**During Manuscript Submission:**
1. Create GitHub release v1.0.0
2. Zenodo auto-creates record (with temporary DOI)
3. Verify and publish on Zenodo
4. Update README.md with Zenodo DOI
5. Reference Zenodo DOI in manuscript's Data Availability section

**After Manuscript Acceptance:**
1. Update .zenodo.json with actual paper DOI
2. Create new GitHub release v1.1.0
3. Zenodo creates new version with updated metadata
4. Use updated version DOI in final paper

---

## Summary: Which Files Matter for Zenodo?

### Critical for Zenodo
- ✓ `.zenodo.json` - **AUTO-FILLS METADATA**
- ✓ `CITATION.cff` - Backup metadata
- ✓ `README.md` - Landing page documentation

### Important for Users
- ✓ `AUTHORS.md` - Attribution
- ✓ `QUICKSTART.md` - Getting started
- ✓ `CHANGELOG.md` - Version history
- ✓ All data and code files

### Not Directly Used by Zenodo (But Good Practice)
- `ACKNOWLEDGMENTS.md` - For credit/funding info
- `plan.md` - Internal planning document
- `.gitignore` - Git-specific, not needed in Zenodo archive

---

## Need Help?

- **Zenodo Documentation:** https://help.zenodo.org/
- **Zenodo Support:** info@zenodo.org
- **GitHub-Zenodo Guide:** https://guides.github.com/activities/citable-code/

---

*Last updated: January 26, 2025*
