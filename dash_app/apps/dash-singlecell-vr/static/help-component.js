AFRAME.registerComponent('help', {
    schema: {
        show: {type: 'boolean', default: false},
        zDepth: {type: 'number', default: -1}
    },
    init: function() {
        this.el.object3D.position.set(0, 0, this.data.zDepth);
        this.width = this.visibleWidthAtZDepth(this.data.zDepth);
        this.height = this.visibleHeightAtZDepth(this.data.zDepth);
        this.bufferWidth = this.width * 3;
        this.bufferHeight = this.height * 3;
        const geometry = new THREE.PlaneBufferGeometry( this.width, this.height );
        const material = new THREE.MeshBasicMaterial( {color: 0x1F2630, side: THREE.DoubleSide} );
        const plane = new THREE.Mesh( geometry, material );
        this.buffer = document.getElementById("helpBuffer");
        this.buffer.object3D.add( plane );
        this.buffer.setAttribute('visible', false);
        this.el.setAttribute('visible', false);
    },
    update: function(oldData) {
        console.log(this.data.show)
        if (this.data.show) {
            document.getElementById("hud").setAttribute('visible', false);
            this.el.setAttribute('visible', true);
            const width = this.visibleWidthAtZDepth(this.data.zDepth);
            const height = this.visibleHeightAtZDepth(this.data.zDepth);
            if (height > width) {
                this.el.setAttribute('width', width);
                this.el.setAttribute('height', height);
                this.buffer.setAttribute('geometry', 'width', width * 3);
                this.buffer.setAttribute('geometry', 'height', height * 3);
                this.buffer.setAttribute('material', 'color', 0x1F2630);
                this.buffer.setAttribute('visible', true);
                this.buffer.object3D.position.set(this.el.getAttribute('position').x, 
                                                this.el.getAttribute('position').y, 
                                                this.el.getAttribute('position').z - 1);
            } else {
                this.el.setAttribute('width', width);
                this.el.setAttribute('height', height);
                this.buffer.setAttribute('visible', false);
            }
            this.el.setAttribute('src', '/assets/tips_1.m4v');
            setTimeout(() => {
                this.el.setAttribute('src', '/assets/tips_2.m4v');
                setTimeout(() => {
                    this.el.setAttribute('src', '/assets/tips_3.m4v');
                    setTimeout(() => {
                        this.dismiss();
                    }, 5000)
                }, 5000);
            }, 5000);
            
        } else {
            this.dismiss();
        }
    },
    // These are duplicated, fix later
    visibleHeightAtZDepth: function( depth ) {
        const camera = this.el.sceneEl.camera;
        // compensate for cameras not positioned at z=0
        const cameraOffset = camera.position.z;
        if ( depth < cameraOffset ) depth -= cameraOffset;
        else depth += cameraOffset;
      
        // vertical fov in radians
        const vFOV = camera.fov * Math.PI / 180;
      
        // Math.abs to ensure the result is always positive
        return 2 * Math.tan( vFOV / 2 ) * Math.abs( depth );
    },
    visibleWidthAtZDepth: function( depth ) {
        const camera = this.el.sceneEl.camera;
        const height = this.visibleHeightAtZDepth( depth, camera );
        let width = height * camera.aspect;
        return width;
    },
    dismiss: function() {
        this.el.setAttribute('visible', false);
        this.buffer.setAttribute('visible', false);
        document.getElementById("hud").setAttribute('visible', true);
    }

});