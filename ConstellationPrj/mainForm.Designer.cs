namespace ConstellationPrj
{
    partial class MainForm
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.Go = new System.Windows.Forms.Button();
            this.trjCompProgressBar = new System.Windows.Forms.ProgressBar();
            this.SuspendLayout();
            // 
            // Go
            // 
            this.Go.Location = new System.Drawing.Point(153, 168);
            this.Go.Name = "Go";
            this.Go.Size = new System.Drawing.Size(75, 23);
            this.Go.TabIndex = 0;
            this.Go.Text = "Start";
            this.Go.UseVisualStyleBackColor = true;
            this.Go.Click += new System.EventHandler(this.Go_Click);
            // 
            // trjCompProgressBar
            // 
            this.trjCompProgressBar.Location = new System.Drawing.Point(93, 92);
            this.trjCompProgressBar.Name = "trjCompProgressBar";
            this.trjCompProgressBar.Size = new System.Drawing.Size(195, 22);
            this.trjCompProgressBar.TabIndex = 1;
            // 
            // MainForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(384, 235);
            this.Controls.Add(this.trjCompProgressBar);
            this.Controls.Add(this.Go);
            this.Name = "MainForm";
            this.Text = "Aircraft trajectory computation";
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.Button Go;
        private System.Windows.Forms.ProgressBar trjCompProgressBar;
    }
}

